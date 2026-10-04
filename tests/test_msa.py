"""Testes offline de MSA real, consenso, conservacao e mapeamento de colunas.

Nenhum teste inventa alinhamento. Fixtures FASTA alinhadas sao entradas de
parser, nao respostas ao vivo da EMBL-EBI.
"""

from __future__ import annotations

import io
import subprocess
from email.message import EmailMessage
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request

import pytest

from modules import alignment, dna_analysis, msa, provenance, scientific_checks

FIXTURES = Path(__file__).parent / "fixtures"
EMAIL = "user@example.com"


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


class _FakeHandle:
    def __init__(self, body: str) -> None:
        self._body = body.encode("utf-8")

    def read(self) -> bytes:
        return self._body

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _members(*sequences: str, molecule: str = "DNA") -> list:
    out = []
    for index, seq in enumerate(sequences, start=1):
        out.append(
            msa.build_collection_member(
                sequence=seq,
                identifier=f"seq_{index}",
                source="user input",
                molecule=molecule,
            )
        )
    return out


def test_pairwise_availability_says_nw_is_not_msa():
    info = alignment.multiple_alignment_availability()
    assert info["available"] is True
    assert info["blast_is_not_msa"] is True
    assert "not MSA" in info["reason"] or "Needleman" in info["reason"]
    blast = alignment.blast_search_availability()
    assert blast["available"] is False


def test_local_tools_absent_are_not_faked(monkeypatch):
    monkeypatch.setattr(msa.shutil, "which", lambda _name: None)
    local = msa.detect_local_tools()
    assert local["mafft"]["available"] is False
    assert local["mafft"]["version"] == ""
    assert "not installed" in local["mafft"]["reason"]
    assert "ebi_clustalo" in msa.available_backend_ids()
    assert "local_mafft" not in msa.available_backend_ids()


def test_parse_fasta_records_keeps_all_entries():
    parsed = dna_analysis.parse_fasta_records(">one\nAAAA\n>two\nCCCC\n")
    assert parsed["record_count"] == 2
    assert parsed["records"][0]["identifier"] == "one"
    assert parsed["records"][1]["sequence"] == "CCCC"
    with pytest.raises(ValueError, match="2 records"):
        dna_analysis.parse_sequence_payload(">one\nAAAA\n>two\nCCCC\n")


def test_mixed_molecules_and_single_sequence_rejected():
    dna = msa.build_collection_member(sequence="ACGTACGT", identifier="d")
    protein = msa.build_collection_member(sequence="MKTFFVAA", identifier="p")
    rna = msa.build_collection_member(sequence="ACGUACGU", identifier="r")
    with pytest.raises(msa.MsaError) as mixed:
        msa.validate_collection([dna, protein])
    assert mixed.value.category == "INVALID_INPUT"
    with pytest.raises(msa.MsaError) as mixed_na:
        msa.validate_collection([dna, rna])
    assert mixed_na.value.category == "INVALID_INPUT"
    with pytest.raises(msa.MsaError) as single:
        msa.validate_collection([dna])
    assert single.value.category == "INVALID_INPUT"


def test_duplicates_are_flagged_not_removed():
    members = _members("ACGTACGT", "ACGTACGT")
    groups = msa.identical_sequence_groups(members)
    assert len(groups) == 1
    assert groups[0]["n"] == 2
    validated = msa.validate_collection(members)
    assert validated["n_sequences"] == 2


def test_resource_limit_sequence_count():
    members = _members("ACGT", "ACGA")
    too_many = members * 20
    with pytest.raises(msa.MsaError) as exc:
        msa.validate_collection(too_many)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_parse_aligned_and_mapping_with_gaps():
    members = _members("ACGT", "ACCGT")
    rows = msa.parse_aligned_fasta(_load("msa_aligned_gaps.fa"))
    validated = msa.validate_aligned_rows(rows, members)
    assert validated["alignment_length"] == 5
    maps = msa.column_coordinate_maps(validated["rows"])
    assert maps[0] == [0, 1, None, 2, 3]
    assert maps[1] == [0, 1, 2, 3, 4]


def test_leading_internal_terminal_gap_mapping():
    rows = [
        {"aligned": "--AC", "ungapped": "AC", "identifier": "a"},
        {"aligned": "A--C", "ungapped": "AC", "identifier": "b"},
        {"aligned": "AC--", "ungapped": "AC", "identifier": "c"},
    ]
    maps = msa.column_coordinate_maps(rows)
    assert maps[0] == [None, None, 0, 1]
    assert maps[1] == [0, None, None, 1]
    assert maps[2] == [0, 1, None, None]


def test_unequal_aligned_rows_are_parsing_error():
    members = _members("ACGT", "ACGTTTTT")
    rows = msa.parse_aligned_fasta(_load("msa_aligned_unequal.fa"))
    with pytest.raises(msa.MsaError) as invalid:
        msa.validate_aligned_rows(rows, members)
    assert invalid.value.category == "PARSING_ERROR"


def test_ungapped_mismatch_is_not_repaired():
    members = _members("AAAA", "CCCC")
    rows = msa.parse_aligned_fasta(_load("msa_aligned_equal.fa"))
    with pytest.raises(msa.MsaError) as exc:
        msa.validate_aligned_rows(rows, members)
    assert exc.value.category == "PARSING_ERROR"


def test_consensus_majority_ties_gaps_and_threshold():
    identical = msa.consensus_from_columns(["AAAA", "AAAA"], "DNA")
    assert identical["sequence"] == "AAAA"
    majority = msa.consensus_from_columns(["AAAA", "AAAT"], "DNA")
    assert majority["sequence"] == "AAAW"
    assert identical["status"] == "COMPUTED"
    assert "not an" in identical["disclaimer"].lower() or "not an experimental" in identical["disclaimer"].lower()
    tie = msa.consensus_from_columns(["A", "T"], "DNA")
    assert tie["sequence"] == "W"
    rna_tie = msa.consensus_from_columns(["A", "U"], "RNA")
    assert rna_tie["sequence"] == "N"
    protein_tie = msa.consensus_from_columns(["A", "V"], "PROTEIN")
    assert protein_tie["sequence"] == "X"
    gapped = msa.consensus_from_columns(["A-", "A-"], "DNA")
    assert gapped["sequence"] == "A-"
    weak = msa.consensus_from_columns(["AAAA", "TTTT"], "DNA", majority_threshold=0.6)
    assert weak["sequence"] == "NNNN"
    with pytest.raises(msa.MsaError):
        msa.consensus_from_columns([], "DNA")
    single = msa.consensus_from_columns(["ACGT"], "DNA")
    assert single["sequence"] == "ACGT"


def test_conservation_identical_variable_and_all_gap():
    ident = msa.conservation_shannon(["AAAA", "AAAA"], "DNA")
    assert ident["scores"][0] == pytest.approx(1.0)
    assert scientific_checks.conservation_score_is_valid(ident["scores"][0])
    variable = msa.conservation_shannon(["AAAA", "TTTT"], "DNA")
    assert variable["scores"][0] == pytest.approx(0.5)
    gapped = msa.conservation_shannon(["A-", "T-"], "DNA")
    assert gapped["gap_fractions"][1] == pytest.approx(1.0)
    assert gapped["scores"][1] != gapped["scores"][1]  # NaN
    assert scientific_checks.conservation_score_is_valid(gapped["scores"][1])
    classes = msa.classify_column_variation(["AA-", "AT-"])
    assert classes[0] == "conserved"
    assert classes[1] == "variable"
    assert classes[2] == "gap"


def test_variation_and_identity_from_msa_not_nw():
    rows = ["AC-GT", "ACCGT"]
    table = msa.variation_vs_reference(rows, 0)
    assert table[1][2] == "insertion"
    ident = msa.identity_matrix_from_msa(["ACGT", "ACGT"])
    assert ident["not_needleman_wunsch"] is True
    assert ident["identity_scale"] == "fraction_0_1"
    assert ident["matrix"][0][1] == pytest.approx(1.0)
    empty_denom = msa.identity_matrix_from_msa(["A-", "-A"])
    assert empty_denom["matrix"][0][1] != empty_denom["matrix"][0][1]


def test_ncbi_record_to_msa_preserves_hash():
    record = {
        "accession": "NM_000000.1",
        "description": "fixture",
        "organism": "Homo sapiens",
        "sequence": "ACGTACGT",
        "sequence_available": True,
        "database": "nucleotide",
        "source": "NCBI Entrez",
        "retrieved_at_utc": "2026-01-01T00:00:00Z",
    }
    member = msa.member_from_ncbi_record(record)
    assert member["sequence"] == "ACGTACGT"
    assert member["hash"] == provenance.sequence_digest("ACGTACGT")
    assert member["accession"] == "NM_000000"
    assert member["version"] == "NM_000000.1"
    assert member["source"] == "NCBI Entrez"


def test_blast_hit_is_not_a_sequence_and_timeout_is_refused():
    hit = {
        "accession": "NM_000000",
        "organism": "Homo sapiens",
        "hseq": "ACGTACGTACGTACGT",
        "qseq": "ACGTACGT",
    }
    with pytest.raises(msa.MsaError) as timed:
        msa.pending_member_from_blast_hit(hit, {"status": "TIMEOUT", "rid": "ABCDEFGH"})
    assert timed.value.category == "INVALID_INPUT"
    pending = msa.pending_member_from_blast_hit(
        hit,
        {
            "status": "READY",
            "rid": "ABCDEFGH",
            "program": "blastn",
            "database": "core_nt",
            "query_hash": "abc",
        },
    )
    assert pending["sequence"] == ""
    assert pending["status"] == "SEQUENCE_UNAVAILABLE"
    assert pending["sequence"] != hit["hseq"]
    record = {
        "accession": "NM_000000.1",
        "sequence": "ACGTACGT",
        "sequence_available": True,
        "database": "nucleotide",
        "source": "NCBI Entrez",
        "retrieved_at_utc": "2026-01-01T00:00:00Z",
        "organism": "Homo sapiens",
    }
    member = msa.complete_blast_member_with_ncbi(pending, record)
    assert member["sequence"] == "ACGTACGT"
    assert member["blast_rid"] == "ABCDEFGH"
    assert member["source"] == "BLAST-derived selection"
    assert member["hash"] == provenance.sequence_digest("ACGTACGT")


def test_ebi_submit_poll_retrieve_uses_real_parser_path():
    aligned = _load("msa_aligned_gaps.fa")
    members = _members("ACGT", "ACCGT")

    def fake(request: Request, timeout=None):
        url = str(request.full_url)
        assert url.startswith(msa.EBI_CLUSTALO_BASE)
        if url.endswith("/run"):
            return _FakeHandle("clustalo-I20260101-000000-0000-12345678-p1m")
        if "/status/" in url:
            return _FakeHandle("FINISHED")
        if "/result/" in url:
            assert "/aln-fasta" in url or url.endswith("/fa") or url.endswith("/out")
            return _FakeHandle(aligned)
        raise AssertionError(url)

    job = msa.submit_msa(members, EMAIL, backend="ebi_clustalo", urlopen_fn=fake)
    assert job["job_id"].startswith("clustalo-")
    job["submitted_monotonic"] = 0.0
    polled = msa.poll_msa(job, EMAIL, urlopen_fn=fake, now=msa.MIN_POLL_INTERVAL_S + 1)
    assert polled["status"] == "COMPLETED"
    result = msa.retrieve_msa(polled, EMAIL, urlopen_fn=fake)
    assert result["alignment_length"] == 5
    assert result["consensus"]["status"] == "COMPUTED"
    assert result["tool"] == "Clustal Omega"
    assert result["input_hashes"] == [m["hash"] for m in members]
    assert result["input_order"] == [m["identifier"] for m in members]
    assert result["parameters"]["order"] == "input"
    assert result["retrieved_at_utc"]
    assert result["tool_version"] == ""
    assert result["status"] == "COMPLETED"
    assert "not Needleman" in result["identity_matrix"]["method"] or result["identity_matrix"]["not_needleman_wunsch"]
    csv_text = msa.export_column_csv(result)
    assert "conservation_shannon" in csv_text
    fasta = msa.export_aligned_fasta(result)
    assert "input_hashes=" in fasta
    detail = msa.column_detail(result, 2)
    assert detail["residues"] == ["-", "C"]
    assert detail["original_positions_0based"][0] is None
    span = msa.spans_for_column(result, 0, 0)
    assert span["sequence"] == "A"
    cached = msa.mark_cached_result(result)
    assert cached["cache_status"] == "cached"
    assert cached["retrieved_at_utc"] == result["retrieved_at_utc"]


def test_ebi_http_429_and_invalid_job_id():
    def boom(request, timeout=None):
        raise HTTPError(
            msa.EBI_CLUSTALO_BASE + "/run",
            429,
            "Too Many Requests",
            hdrs=EmailMessage(),
            fp=io.BytesIO(b""),
        )

    members = _members("ACGTACGT", "ACGTACGA")
    with pytest.raises(msa.MsaError) as exc:
        msa.submit_msa(members, EMAIL, urlopen_fn=boom)
    assert exc.value.category == "RATE_LIMITED"
    with pytest.raises(msa.MsaError) as job:
        msa.poll_msa(
            {
                "backend": "ebi_clustalo",
                "job_id": "https://evil.example/msa",
                "submitted_monotonic": 0.0,
            },
            EMAIL,
            now=10_000,
        )
    assert job.value.category == "INVALID_INPUT"


def test_poll_timeout_is_not_unavailable():
    with pytest.raises(msa.MsaError) as exc:
        msa.poll_msa(
            {
                "backend": "ebi_clustalo",
                "job_id": "clustalo-I20260101-000000-0000-12345678-p1m",
                "submitted_monotonic": 0.0,
                "last_poll_monotonic": 0.0,
            },
            EMAIL,
            now=msa.MAX_JOB_WAIT_S + 1,
        )
    assert exc.value.category == "TIMEOUT"


def test_tool_not_installed_local_backend(monkeypatch):
    monkeypatch.setattr(msa.shutil, "which", lambda _name: None)
    members = _members("ACGTACGT", "ACGTACGA")
    with pytest.raises(msa.MsaError) as exc:
        msa.submit_msa(members, EMAIL, backend="local_mafft")
    assert exc.value.category == "TOOL_NOT_INSTALLED"


def test_malformed_alignment_empty_is_parsing_error():
    with pytest.raises(msa.MsaError) as exc:
        msa.parse_aligned_fasta("   ")
    assert exc.value.category == "PARSING_ERROR"


def test_cache_key_includes_order_and_backend():
    hashes = ["aaa", "bbb"]
    a = msa.cache_key(hashes=hashes, backend="ebi_clustalo", tool_version="", parameters={"order": "input"})
    b = msa.cache_key(hashes=["bbb", "aaa"], backend="ebi_clustalo", tool_version="", parameters={"order": "input"})
    c = msa.cache_key(hashes=hashes, backend="local_mafft", tool_version="", parameters={"order": "input"})
    d = msa.cache_key(hashes=hashes, backend="ebi_clustalo", tool_version="1.2.4", parameters={"order": "input"})
    assert a != b
    assert a != c
    assert a != d


def test_scientific_checks_for_msa_scores():
    assert scientific_checks.msa_row_length_matches("AC-GT", 5)
    assert scientific_checks.msa_row_length_matches("AC-GT", 4) is False
    assert scientific_checks.conservation_score_is_valid(0.0)
    assert scientific_checks.conservation_score_is_valid(1.0)
    assert scientific_checks.conservation_score_is_valid(float("nan"))
    assert scientific_checks.conservation_score_is_valid(1.5) is False
    assert scientific_checks.conservation_score_is_valid(-0.1) is False


def test_output_size_limit_is_resource_limit():
    huge = "x" * (msa.MAX_OUTPUT_BYTES + 1)
    with pytest.raises(msa.MsaError) as exc:
        msa.parse_aligned_fasta(huge)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_total_residues_limit(monkeypatch):
    monkeypatch.setattr(msa, "MAX_RESIDUES_TOTAL", 10)
    members = _members("ACGTACGT", "ACGTACGA")
    with pytest.raises(msa.MsaError) as exc:
        msa.validate_collection(members)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_local_timeout_is_timeout_not_unavailable(monkeypatch, tmp_path):
    fake_bin = tmp_path / "mafft"
    fake_bin.write_text("", encoding="utf-8")
    monkeypatch.setattr(msa, "_probe_local_version", lambda *_args, **_kwargs: "v1")

    def boom(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd="mafft", timeout=1)

    monkeypatch.setattr(msa.subprocess, "run", boom)
    with pytest.raises(msa.MsaError) as exc:
        msa._execute_local(str(fake_bin), "mafft", ">a\nACGT\n>b\nACGA\n", "DNA")
    assert exc.value.category == "TIMEOUT"


def test_muscle_args_follow_detected_version():
    v3 = msa._local_args("/bin/muscle", "muscle", "in.fa", "out.fa", "DNA", version="MUSCLE 3.8.31")
    assert v3[1:5] == ["-in", "in.fa", "-out", "out.fa"]
    v5 = msa._local_args("/bin/muscle", "muscle", "in.fa", "out.fa", "DNA", version="MUSCLE v5.1")
    assert v5[1:5] == ["-align", "in.fa", "-output", "out.fa"]
    clustal = msa._local_args("/bin/clustalo", "clustalo", "in.fa", "out.fa", "DNA")
    assert clustal[0] == "/bin/clustalo"
    assert all(not arg.startswith(";") for arg in clustal)


def test_ebi_version_is_not_hardcoded():
    info = msa.tool_availability()
    ebi = next(item for item in info["backends"] if item["id"] == "ebi_clustalo")
    assert ebi["version"] == ""
    assert ebi["available"] is True


def test_allowlisted_basename_accepts_windows_exe():
    allowed = ("mafft", "mafft.exe")
    assert msa._normalized_tool_basename(r"C:\bin\mafft.exe") == "mafft"
    assert "mafft.exe" in msa._allowlisted_basenames(allowed)
    assert msa._normalized_tool_basename(r"C:\bin\evil.exe") == "evil"


def test_ebi_result_falls_back_from_aln_fasta_404():
    aligned = _load("msa_aligned_equal.fa")
    members = _members("ACGT", "ACGT")
    calls: list[str] = []

    def fake(request: Request, timeout=None):
        url = str(request.full_url)
        calls.append(url)
        if url.endswith("/run"):
            return _FakeHandle("clustalo-I20260101-000000-0000-12345678-p1m")
        if "/status/" in url:
            return _FakeHandle("FINISHED")
        if url.endswith("/aln-fasta"):
            raise HTTPError(
                url,
                404,
                "Not Found",
                hdrs=EmailMessage(),
                fp=io.BytesIO(b""),
            )
        if "/result/" in url:
            return _FakeHandle(aligned)
        raise AssertionError(url)

    job = msa.submit_msa(members, EMAIL, backend="ebi_clustalo", urlopen_fn=fake)
    job["submitted_monotonic"] = 0.0
    polled = msa.poll_msa(job, EMAIL, urlopen_fn=fake, now=msa.MIN_POLL_INTERVAL_S + 1)
    result = msa.retrieve_msa(polled, EMAIL, urlopen_fn=fake)
    assert result["alignment_length"] == 4
    assert any(url.endswith("/aln-fasta") for url in calls)
    assert any(url.endswith("/fa") or url.endswith("/out") for url in calls)


def test_msa_envelope_includes_alignment_hash():
    members = _members("ACGT", "ACGT")
    result = msa.build_msa_result(
        raw_alignment=_load("msa_aligned_equal.fa"),
        members=members,
        tool="fixture",
        tool_version="test",
        method="prealigned fixture",
        parameters={},
    )
    assert result["alignment_hash"]
    assert len(result["alignment_hash"]) == 64
    assert "phylogenetic tree" in result["disclaimer"].lower()
