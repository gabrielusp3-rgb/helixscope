"""Fase 15: OpenCL gates, completeness, FULL GENOME lock, annotation, worker."""

from __future__ import annotations

import os
from types import SimpleNamespace

from modules import (
    crispr,
    crispr_assemblies,
    crispr_casoffinder,
    crispr_offtarget,
    engine_validation,
    genome_jobs,
    genome_region,
    genome_store,
    provenance,
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


def test_parse_cas_offinder_header_with_spaces():
    text = (
        "GAGTCCGAGCAGAAGAAGAANGG\tHS_TEST_1 HelixScope TEST REFERENCE contig_1; "
        "not GRCh38\t40\tGAGTCCGAGCAGAAGAAGAAGGG\t+\t0\n"
    )
    parsed = crispr_casoffinder.parse_cas_offinder_output(text)
    assert parsed["status"] == "OK"
    assert parsed["rows"][0]["chromosome"] == "HS_TEST_1"
    assert parsed["rows"][0]["position_reported"] == 40
    assert parsed["rows"][0]["strand"] == "+"
    assert parsed["rows"][0]["mismatches_reported"] == 0
    policy = crispr_casoffinder.cas_offinder_version_policy()
    assert policy["production_version"] == "2.4.1"
    assert policy["v3_production_ready"] is False
    assert policy["native_bulges"] is False
    assert policy["bulge_status"] == "UNAVAILABLE"


def test_bulges_are_unavailable_not_cas3(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        dna_bulge=1,
        rna_bulge=0,
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(""),
    )
    assert envelope["status"] == "UNAVAILABLE"
    assert envelope["genome_wide"] is False
    assert "3" in envelope["reason"] or "bulge" in envelope["reason"].lower()


def test_full_genome_blocked_without_ready_or_live(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    engine_validation.reset_live_validations()
    gate = genome_jobs.full_genome_search_preflight("GRCh38.p14")
    assert gate["allowed"] is False
    assert "READY" in gate["reason"] or "Download" in gate["reason"]
    fasta = tmp_path / "refs" / "GRCh38.p14" / "genome.fa"
    fasta.parent.mkdir(parents=True)
    fasta.write_bytes(b">NC_TEST\n" + b"ACGT" * 40 + b"\n")
    genome_store.finalize_ready_from_fasta(
        "GRCh38.p14",
        str(fasta),
        kind="public_assembly",
        not_a_public_assembly=False,
        official_compressed_md5="c30471567037b2b2389d43c908c653e1",
        compressed_md5_observed="c30471567037b2b2389d43c908c653e1",
    )
    gate2 = genome_jobs.full_genome_search_preflight("GRCh38.p14")
    assert gate2["allowed"] is False
    reason = str(gate2["reason"])
    assert (
        "LIVE_VALIDATED" in reason
        or "OpenCL" in reason
        or "not installed" in reason.lower()
    )
    try:
        genome_jobs.submit_casoffinder_job(
            guide=_emx1_guide(),
            assembly_id="GRCh38.p14",
            spawn=False,
        )
        assert False, "public submit must be blocked without LIVE_VALIDATED"
    except genome_jobs.GenomeJobError as exc:
        assert exc.category == "UNAVAILABLE"


def test_test_reference_complete_is_never_genome_wide(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    line = f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        include_nag=False,
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(line),
    )
    assert envelope["status"] == crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE
    assert envelope["genome_wide"] is False
    assert envelope["search_completeness"]["complete"] is True
    assert envelope["engine_input_manifest"]["sha256_match"] is True
    assert envelope["guide_specificity"]["mit_sguide"]["score"] is not None
    assert envelope["guide_specificity"]["mit_sguide"]["genome_wide"] is False


def test_hit_cap_is_not_completed_full_reference(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    line = f"{EMX1_SITE} HS_TEST_1 40 {EMX1_SITE} + 0\n"
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        include_nag=False,
        max_hits=1,
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(line + line),
    )
    assert envelope["status"] == crispr_offtarget.JOB_RESOURCE_LIMIT
    assert envelope["genome_wide"] is False
    assert envelope["search_completeness"]["complete"] is False
    assert envelope["guide_specificity"]["mit_sguide"]["score"] is None
    assert envelope["guide_specificity"]["cfd_sguide"]["score"] is None


def test_sha256_mismatch_refuses_search(tmp_path, monkeypatch):
    _isolate_store(tmp_path, monkeypatch)
    genome_store.install_test_reference(fixture_path=FIXTURE)
    record = genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256="00" * 32,
        include_nag=False,
        which_info=_which_ok(tmp_path),
        run_fn=_fake_casoffinder_run(""),
    )
    assert envelope["status"] == crispr_offtarget.JOB_FAILED
    assert envelope["genome_wide"] is False


def test_completeness_helpers():
    manifest = {
        "sha256_match": True,
        "n_contigs": 2,
    }
    ok = crispr_casoffinder.search_completeness_record(
        manifest=manifest,
        truncated=False,
        timed_out=False,
        cancelled=False,
        tool_error="",
        exit_ok=True,
    )
    assert ok["complete"] is True
    capped = crispr_casoffinder.search_completeness_record(
        manifest=manifest,
        truncated=True,
        timed_out=False,
        cancelled=False,
        tool_error="",
        exit_ok=True,
    )
    assert capped["complete"] is False
    timed = crispr_casoffinder.search_completeness_record(
        manifest=manifest,
        truncated=False,
        timed_out=True,
        cancelled=False,
        tool_error="",
        exit_ok=False,
    )
    assert timed["complete"] is False


def test_region_view_is_local_not_full_genome():
    hit = {
        "start_0based": 40,
        "end_0based": 60,
        "pam_start_0based": 60,
        "pam_end_0based": 63,
        "strand": "+",
        "target_sequence": EMX1_SITE,
        "PAM": "GGG",
    }
    view = genome_region.local_region_view(
        contig="HS_TEST_1",
        contig_length=200,
        hit=hit,
        features=[
            {
                "seqid": "HS_TEST_1",
                "type": "gene",
                "start_0based": 40,
                "end_0based": 80,
                "strand": "+",
                "gene_id": "HSG1",
                "gene_name": "TESTGENE",
            }
        ],
        flank_nt=10,
    )
    assert view["full_genome_dom"] is False
    assert view["axis"]["end_0based"] - view["axis"]["start_0based"] == 40
    shifted = genome_region.shift_window(view, delta_nt=-5, contig_length=200)
    assert shifted["start_0based"] == view["axis"]["start_0based"] - 5


def test_job_log_sanitizes_home(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIXSCOPE_JOBS_DIR", str(tmp_path / "jobs"))
    home = os.path.expanduser("~")
    cleaned = genome_jobs.sanitize_job_log(f"path {home}\\secret\\file.txt NCBI_KEY")
    assert home not in cleaned
    assert "<home>" in cleaned


def test_reference_scope_not_silent_primary_only():
    human = crispr_assemblies.catalog_by_id("GRCh38.p14")
    assert "primary" in human["reference_scope"].lower()
    assert "silently" in human["reference_scope"].lower()
    assert human["accession"] == "GCF_000001405.40"


def test_live_validation_persists_across_memory_reset(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIXSCOPE_JOBS_DIR", str(tmp_path / "jobs"))
    engine_validation.reset_live_validations()
    engine_validation.record_live_validation("Cas-OFFinder", ok=True, version="2.4.1-test")
    engine_validation.reset_live_validations()
    loaded = engine_validation.live_record("Cas-OFFinder")
    assert loaded is not None
    assert loaded["status"] == engine_validation.STATUS_LIVE_VALIDATED
    assert provenance.HELIXSCOPE_VERSION.startswith("0.")


def test_full_genome_blocked_when_memory_insufficient(tmp_path, monkeypatch):
    from modules import resource_admission

    _isolate_store(tmp_path, monkeypatch)
    engine_validation.reset_live_validations()
    engine_validation.record_live_validation(
        "Cas-OFFinder", ok=True, version="2.4.1-test"
    )
    from modules import crispr_casoffinder

    monkeypatch.setattr(
        crispr_casoffinder,
        "detect_cas_offinder",
        lambda: {"available": True, "status": "installed", "version": "2.4.1-test"},
    )
    fasta = tmp_path / "refs" / "GRCh38.p14" / "genome.fa"
    fasta.parent.mkdir(parents=True)
    fasta.write_bytes(b">NC_TEST\n" + b"ACGT" * 40 + b"\n")
    genome_store.finalize_ready_from_fasta(
        "GRCh38.p14",
        str(fasta),
        kind="public_assembly",
        not_a_public_assembly=False,
        official_compressed_md5="c30471567037b2b2389d43c908c653e1",
        compressed_md5_observed="c30471567037b2b2389d43c908c653e1",
    )

    def refuse(**kwargs):
        return {
            "decision": resource_admission.DECISION_RESOURCE_LIMIT,
            "reason": "Injected: not enough RAM for this reference.",
            "blocker": "memory",
        }

    monkeypatch.setattr(resource_admission, "admit_reference_search", refuse)
    gate = genome_jobs.full_genome_search_preflight("GRCh38.p14")
    assert gate["allowed"] is False
    assert gate.get("resource_status") == resource_admission.DECISION_RESOURCE_LIMIT
    assert "RAM" in gate["reason"] or "memory" in gate["reason"].lower()
