"""Fase 14: referencias verificadas, FAI, jobs, Cas-OFFinder, gates de score."""

from __future__ import annotations

import gzip
import hashlib
import os
import zipfile
from types import SimpleNamespace

import pytest

from modules import (
    crispr,
    crispr_assemblies,
    crispr_casoffinder,
    crispr_offtarget,
    crispr_specificity,
    genome_download,
    genome_fasta,
    genome_jobs,
    genome_store,
    genome_worker,
)


EMX1_SPACER = "GAGTCCGAGCAGAAGAAGAA"
EMX1_SITE = crispr.EXAMPLE_EMX1_SPCAS9
FIXTURE = os.path.join(
    os.path.dirname(__file__), "fixtures", "helixscope_test_reference.fa"
)


def _emx1_guide() -> dict:
    guides = crispr.find_guides(EMX1_SITE, "SpCas9")
    return next(g for g in guides if g["guide_sequence"] == EMX1_SPACER)


def _isolate_store(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("HELIXSCOPE_REFERENCE_DIR", str(tmp_path / "refs"))
    monkeypatch.setenv("HELIXSCOPE_JOBS_DIR", str(tmp_path / "jobs"))


def _fake_casoffinder_run(output_text: str):
    def run_fn(args, **kwargs):
        out = args[3]
        with open(out, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(output_text)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    return run_fn


def _which_ok(tmp_path) -> dict:
    exe = tmp_path / "cas-offinder.exe"
    exe.write_text("not-a-real-binary", encoding="utf-8")
    return {
        "available": True,
        "path": str(exe),
        "version": "2.4.1-test",
        "reason": "test double",
        "report_path": "cas-offinder.exe",
    }


def test_catalog_has_official_ncbi_md5_not_placeholder():
    human = crispr_assemblies.catalog_by_id("GRCh38.p14")
    t2t = crispr_assemblies.catalog_by_id("T2T-CHM13v2.0")
    mouse = crispr_assemblies.catalog_by_id("GRCm39")
    assert human["official_compressed_md5"] == "c30471567037b2b2389d43c908c653e1"
    assert t2t["official_compressed_md5"] == "9e6bf6b586bc8954208d1cc1d5f2fc99"
    assert mouse["official_compressed_md5"] == "c0b0c4c3f54d2b480efe68a18bf7e42b"
    assert human["checksum"] == ""
    assert "placeholder" not in human["official_compressed_md5"]
    test_ref = crispr_assemblies.test_reference_descriptor()
    assert test_ref["kind"] == "synthetic_test_reference"
    assert test_ref["not_a_public_assembly"] is True


def test_gtf_grch37_rejected_on_grch38():
    assert crispr_assemblies.gtf_matches_assembly("GRCh38.p14", "GRCh38.p14") is True
    assert crispr_assemblies.gtf_matches_assembly("GRCh38.p14", "hg19") is False
    assert crispr_assemblies.gtf_matches_assembly("GRCh38.p14", "GRCh37") is False
    assert crispr_assemblies.gtf_matches_assembly("GRCh38.p14", "") is False


def test_filename_hg38_is_not_verified_assembly():
    lookup = crispr_assemblies.lookup_declaration("hg38")
    assert lookup["verified_assembly"] == ""
    assert lookup["filename_trusted_as_assembly"] is False
    assert lookup["checksum"] == ""


def test_fai_random_access_coordinates_and_rc(tmp_path):
    fasta = tmp_path / "mini.fa"
    fasta.write_bytes(b">NC_TEST_KEEP\nATGCATGCAT\n>HS_TEST_2\nGGGGCCCC\n")
    indexed = genome_fasta.build_fai_for_fasta(str(fasta))
    assert indexed["n_contigs"] == 2
    assert indexed["reference_hash"] == indexed["file_sha256"]
    fai = genome_fasta.load_fai(indexed["fai_path"])
    assert "NC_TEST_KEEP" in fai["contigs"]
    assert "chr1" not in fai["contigs"]
    row = fai["contigs"]["NC_TEST_KEEP"]
    assert genome_fasta.fetch_sequence(str(fasta), row, 0, 10) == "ATGCATGCAT"
    assert genome_fasta.fetch_sequence(str(fasta), row, 0, 1) == "A"
    assert genome_fasta.fetch_sequence(str(fasta), row, 9, 10) == "T"
    assert genome_fasta.fetch_sequence(str(fasta), row, 0, 4, strand="-") == "GCAT"
    disp = genome_fasta.display_interval_1based_inclusive(0, 10)
    assert disp["start_1based"] == 1
    assert disp["end_1based"] == 10
    empty = genome_fasta.display_interval_1based_inclusive(4, 4)
    assert empty["empty"] is True
    with pytest.raises(genome_fasta.GenomeFastaError):
        genome_fasta.fetch_sequence(str(fasta), row, 0, 11)
    with pytest.raises(genome_fasta.GenomeFastaError):
        genome_fasta.fetch_sequence(str(fasta), row, -1, 2)
    meta = genome_fasta.read_fai_provenance(indexed["fai_path"])
    assert meta["reference_sha256"] == indexed["file_sha256"]
    assert genome_fasta.fai_is_current(
        str(fasta), indexed["fai_path"], indexed["file_sha256"]
    )


def test_fai_stale_when_fasta_size_changes(tmp_path):
    fasta = tmp_path / "mini.fa"
    fasta.write_bytes(b">c1\nACGTACGTAC\n")
    indexed = genome_fasta.build_fai_for_fasta(str(fasta))
    fasta.write_bytes(b">c1\nACGTACGTACNNNN\n")
    assert genome_fasta.fai_is_current(
        str(fasta), indexed["fai_path"], indexed["file_sha256"]
    ) is False


def test_checksum_valid_invalid_missing(tmp_path):
    path = tmp_path / "blob.bin"
    path.write_bytes(b"helixscope-md5-vector")
    digest = genome_fasta.file_md5(str(path))
    assert digest == hashlib.md5(b"helixscope-md5-vector").hexdigest()
    assert digest != hashlib.md5(b"corrupt").hexdigest()
    with pytest.raises(genome_fasta.GenomeFastaError):
        genome_fasta.file_md5(str(tmp_path / "missing.bin"))


def test_partial_download_never_ready(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.mark_partial("GRCh38.p14", reason="incomplete .partial file")
    assert genome_store.ready_record("GRCh38.p14") is None
    manifest = genome_store.read_manifest("GRCh38.p14")
    assert manifest["status"] == genome_store.STATUS_PARTIAL


def test_corrupt_checksum_never_ready(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    fasta = tmp_path / "refs" / "GRCh38.p14" / "genome.fa"
    fasta.parent.mkdir(parents=True)
    fasta.write_bytes(b">c1\nACGT\n")
    with pytest.raises(genome_store.GenomeStoreError):
        genome_store.finalize_ready_from_fasta(
            "GRCh38.p14",
            str(fasta),
            kind="public_assembly",
            not_a_public_assembly=False,
            official_compressed_md5="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            compressed_md5_observed="bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        )
    assert genome_store.read_manifest("GRCh38.p14")["status"] == genome_store.STATUS_CORRUPT
    assert genome_store.ready_record("GRCh38.p14") is None


def test_public_ready_refuses_empty_checksum(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    fasta = tmp_path / "refs" / "GRCh38.p14" / "genome.fa"
    fasta.parent.mkdir(parents=True)
    fasta.write_bytes(b">c1\nACGT\n")
    with pytest.raises(genome_store.GenomeStoreError):
        genome_store.finalize_ready_from_fasta(
            "GRCh38.p14",
            str(fasta),
            kind="public_assembly",
            not_a_public_assembly=False,
            official_compressed_md5="",
            compressed_md5_observed="",
        )


def test_install_test_reference_ready_not_genome_wide(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    manifest = genome_store.install_test_reference(fixture_path=FIXTURE)
    assert manifest["status"] == genome_store.STATUS_READY
    assert manifest["not_a_public_assembly"] is True
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    assert record is not None
    fai = genome_fasta.load_fai(record["fai_path"])
    assert "HS_TEST_1" in fai["contigs"]
    row = fai["contigs"]["HS_TEST_1"]
    site = genome_fasta.fetch_sequence(record["fasta_path"], row, 40, 63)
    assert site == EMX1_SITE
    assert genome_store.any_public_assembly_ready() is False


def test_path_jail_rejects_arbitrary_fasta(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    outside = tmp_path / "evil.fa"
    outside.write_bytes(b">x\nACGT\n")
    assert genome_store.path_is_in_store(str(outside)) is False
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=str(outside),
        fai_path=str(outside) + ".fai",
        reference_sha256="00",
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(""),
    )
    assert envelope["status"] == crispr_offtarget.JOB_INVALID
    assert envelope["genome_wide"] is False


def test_url_allowlist_and_http_rejected():
    assert genome_download.url_is_allowed(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/x.gz"
    )
    assert not genome_download.url_is_allowed(
        "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/x.gz"
    )
    assert not genome_download.url_is_allowed("https://evil.example/genomes/all/x.gz")
    assert not genome_download.url_is_allowed(
        "https://github.com/evil/malware/releases/download/1/x.zip"
    )


def test_zip_path_traversal_rejected(tmp_path):
    archive = tmp_path / "evil.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../evil.exe", b"nope")
    with pytest.raises(genome_download.GenomeDownloadError):
        genome_download.extract_zip_safe(str(archive), str(tmp_path / "out"))


def test_gzip_bomb_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(genome_download, "MAX_UNCOMPRESSED_GENOME_BYTES", 64)
    gz_path = tmp_path / "bomb.gz"
    with gzip.open(gz_path, "wb") as handle:
        handle.write(b"A" * 200)
    with pytest.raises(genome_download.GenomeDownloadError) as exc:
        genome_download.decompress_gzip_to_fasta(str(gz_path), str(tmp_path / "out.fa"))
    assert exc.value.category == "RESOURCE_LIMIT"
    assert not (tmp_path / "out.fa").exists()


def test_download_md5_mismatch_keeps_partial(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)

    class _Resp:
        def read(self, n):
            if getattr(self, "_done", False):
                return b""
            self._done = True
            return b"not-the-genome"

        def close(self):
            return None

    def opener(request, timeout=0):
        return _Resp()

    partial = tmp_path / "refs" / "x.partial"
    partial.parent.mkdir(parents=True)
    with pytest.raises(genome_download.GenomeDownloadError) as exc:
        genome_download.download_to_partial(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/x.gz",
            str(partial),
            expected_md5="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            max_bytes=1024,
            urlopen_fn=opener,
        )
    assert exc.value.category == "CORRUPT"


def test_casoffinder_absent_is_unavailable(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info={"available": False, "reason": "not installed", "path": ""},
    )
    assert envelope["status"] == "UNAVAILABLE"
    assert envelope["genome_wide"] is False
    assert envelope.get("verified_hit_count") is None


def test_casoffinder_malformed_output(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run("not enough fields"),
    )
    assert envelope["status"] == "PARSING_ERROR"
    assert envelope["genome_wide"] is False


def test_casoffinder_test_reference_not_genome_wide(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    line = (
        f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
    )
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(line),
        include_nag=False,
    )
    assert envelope["status"] == crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE
    assert envelope["genome_wide"] is False
    assert envelope["verified_hit_count"] >= 1
    hit = envelope["hits"][0]
    assert hit["chromosome_or_contig"] == "HS_TEST_1"
    assert hit["target_sequence"] == EMX1_SPACER
    assert hit["feature_annotation"] == "UNKNOWN"
    scores = envelope["guide_specificity"]
    assert scores["mit_sguide"]["genome_wide"] is False


def test_hit_cap_is_resource_limit_not_genome_wide(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    line = (
        f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
        f"{EMX1_SITE} HS_TEST_1 103 {EMX1_SITE} + 0\n"
    )
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(line),
        include_nag=False,
        max_hits=1,
    )
    assert envelope["status"] == crispr_offtarget.JOB_RESOURCE_LIMIT
    assert envelope["truncated"] is True
    assert envelope["genome_wide"] is False
    assert envelope["guide_specificity"]["mit_sguide"]["score"] is None


def test_timeout_and_cancel_have_no_aggregate(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)

    def timeout_fn(args, **kwargs):
        import subprocess

        raise subprocess.TimeoutExpired(cmd=args, timeout=1)

    timed = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        run_fn=timeout_fn,
        include_nag=False,
    )
    assert timed["status"] == crispr_offtarget.JOB_TIMEOUT
    assert timed["guide_specificity"]["cfd_sguide"]["score"] is None
    cancelled = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        cancel_check=lambda: True,
    )
    assert cancelled["status"] == crispr_offtarget.JOB_CANCELLED
    assert cancelled["genome_wide"] is False


def test_malicious_guide_rejected(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    evil = dict(_emx1_guide())
    evil["guide_sequence"] = "GAGT; rm -rf /"
    envelope = crispr_casoffinder.search_indexed_fasta(
        evil,
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(""),
    )
    assert envelope["status"] == crispr_offtarget.JOB_INVALID


def test_score_gate_partial_never_zero():
    search = {
        "status": "TIMEOUT",
        "truncated": True,
        "genome_wide": True,
        "hits": [{"site_role": "off_target", "score": 50, "cfd_score": 0.2}],
        "scope": "full_genome",
    }
    scores = crispr_specificity.scores_from_search(search)
    assert scores["mit_sguide"]["score"] is None
    assert scores["cfd_sguide"]["score"] is None
    assert scores["genome_wide"] is False


def test_completed_full_reference_may_be_genome_wide():
    guide = _emx1_guide()
    fasta = ">HS_TEST_1\n" + EMX1_SITE + "\n"
    from modules import crispr_reference

    reference = crispr_reference.load_reference_from_text(fasta, source="test")
    local = crispr_offtarget.search_off_targets_in_reference(guide, reference)
    envelope = crispr_offtarget.build_search_envelope(
        status=crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE,
        reason="test full",
        guide=guide,
        reference={**reference, "verified_assembly": "GRCh38.p14"},
        hits=local["hits"],
        max_mm=3,
        include_perfect=True,
        truncated=False,
        hit_cap=200,
        time_cap=20,
        target_hash="",
        elapsed_ms=1,
        search_ms=1,
        index_hash="abc",
        index_n_sites=1,
        index_build_ms=1,
        algorithm="cas_offinder_spcas9_ngg_nag_no_bulge",
        genome_wide=True,
    )
    assert envelope["genome_wide"] is True
    assert envelope["guide_specificity"]["mit_sguide"]["genome_wide"] is True
    truncated = crispr_offtarget.build_search_envelope(
        status=crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE,
        reason="capped",
        guide=guide,
        reference=reference,
        hits=local["hits"],
        max_mm=3,
        include_perfect=True,
        truncated=True,
        hit_cap=1,
        time_cap=20,
        target_hash="",
        elapsed_ms=1,
        search_ms=1,
        index_hash="abc",
        index_n_sites=1,
        index_build_ms=1,
        algorithm="cas_offinder_spcas9_ngg_nag_no_bulge",
        genome_wide=True,
    )
    assert truncated["genome_wide"] is False
    assert truncated["guide_specificity"]["mit_sguide"]["score"] is None


def test_cache_identity_is_assembly_and_engine_sensitive():
    a = genome_jobs.cache_identity(
        guide_hash="aa",
        assembly_id="GRCh38.p14",
        reference_sha256="111",
        engine="cas_offinder",
        engine_version="2.4.1",
        max_mismatches=3,
        dna_bulge=0,
        rna_bulge=0,
        include_nag=True,
    )
    b = genome_jobs.cache_identity(
        guide_hash="aa",
        assembly_id="GRCh38.p14",
        reference_sha256="222",
        engine="cas_offinder",
        engine_version="2.4.1",
        max_mismatches=3,
        dna_bulge=0,
        rna_bulge=0,
        include_nag=True,
    )
    c = genome_jobs.cache_identity(
        guide_hash="aa",
        assembly_id="GRCh38.p14",
        reference_sha256="111",
        engine="cas_offinder",
        engine_version="2.4.2",
        max_mismatches=3,
        dna_bulge=0,
        rna_bulge=0,
        include_nag=True,
    )
    assert a != b
    assert a != c


def test_worker_complete_cancel_and_dead_pid(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    status = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
    )
    job_id = status["job_id"]
    assert len(job_id) == 32
    line = f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
    original = crispr_casoffinder.search_indexed_fasta

    def wrapped(guide, **kwargs):
        kwargs["run_fn"] = _fake_casoffinder_run(line)
        kwargs["which_info"] = _which_ok(tmp_path)
        kwargs["include_nag"] = False
        return original(guide, **kwargs)

    monkeypatch.setattr(crispr_casoffinder, "search_indexed_fasta", wrapped)
    done = genome_worker.run_casoffinder_job(job_id)
    assert done["status"] == crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE
    result = genome_jobs.load_result(job_id)
    assert result["genome_wide"] is False

    other = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
        max_mismatches=2,
    )
    other_id = other["job_id"]
    os.makedirs(genome_jobs.job_dir(other_id), exist_ok=True)
    open(os.path.join(genome_jobs.job_dir(other_id), "cancel"), "w", encoding="utf-8").write("1")
    cancelled = genome_worker.run_casoffinder_job(other_id)
    assert cancelled["status"] == genome_jobs.JOB_CANCELLED

    dead = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
        max_mismatches=1,
    )
    genome_jobs.write_status(
        dead["job_id"],
        {"status": genome_jobs.JOB_RUNNING, "pid": 999999, "reason": "gone"},
    )
    recovered = genome_jobs.load_status(dead["job_id"])
    assert recovered["status"] == genome_jobs.JOB_INTERRUPTED


def test_job_int_keeps_zero_mismatches():
    assert genome_worker.job_int({"max_mismatches": 0}, "max_mismatches", 3) == 0
    assert genome_worker.job_int({}, "max_mismatches", 3) == 3
    assert genome_worker.job_int({"max_mismatches": ""}, "max_mismatches", 3) == 3


def test_worker_passes_zero_max_mismatches(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    status = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
        max_mismatches=0,
    )
    captured: dict = {}
    original = crispr_casoffinder.search_indexed_fasta

    def wrapped(guide, **kwargs):
        captured["max_mismatches"] = kwargs.get("max_mismatches")
        kwargs["run_fn"] = _fake_casoffinder_run(
            f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
        )
        kwargs["which_info"] = _which_ok(tmp_path)
        kwargs["include_nag"] = False
        return original(guide, **kwargs)

    monkeypatch.setattr(crispr_casoffinder, "search_indexed_fasta", wrapped)
    genome_worker.run_casoffinder_job(status["job_id"])
    assert captured["max_mismatches"] == 0


def test_unsafe_assembly_id_rejected(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    with pytest.raises(genome_store.GenomeStoreError):
        genome_store.assembly_dir("..\\Windows")
    with pytest.raises(genome_jobs.GenomeJobError):
        genome_jobs.job_dir("not-hex")


def test_job_cache_reuse(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    first = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
    )
    job_id = first["job_id"]
    envelope = {
        "status": crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE,
        "truncated": False,
        "genome_wide": False,
        "reason": "cached test",
    }
    genome_jobs._atomic_write(
        os.path.join(genome_jobs.job_dir(job_id), "result.json"),
        envelope,
    )
    genome_jobs.write_status(
        job_id,
        {"status": crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE, "reason": "done"},
    )
    again = genome_jobs.submit_casoffinder_job(
        guide=_emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        spawn=False,
        include_nag=False,
    )
    assert again.get("cache_hit") is True
    assert again.get("job_id") == job_id


def test_no_fake_genome_wide_strings_in_new_modules():
    root = os.path.join(os.path.dirname(__file__), "..", "modules")
    banned = ("fake genome", "dummy assembly", "invented coordinate")
    for name in (
        "genome_store.py",
        "genome_fasta.py",
        "genome_jobs.py",
        "genome_download.py",
        "crispr_casoffinder.py",
    ):
        text = open(os.path.join(root, name), encoding="utf-8").read().lower()
        for token in banned:
            assert token not in text
