"""Live Cas-OFFinder on the HelixScope TEST REFERENCE. Skip if absent.

A READY public assembly live run is gated by HELIXSCOPE_LIVE_GENOME=1 and
is never started from this file automatically (no GRCh38 download).
"""

from __future__ import annotations

import os

import pytest

from modules import (
    crispr,
    crispr_assemblies,
    crispr_casoffinder,
    crispr_offtarget,
    genome_store,
)

EMX1_SITE = crispr.EXAMPLE_EMX1_SPCAS9
FIXTURE = os.path.join(
    os.path.dirname(__file__), "fixtures", "helixscope_test_reference.fa"
)


def _emx1_guide() -> dict:
    guides = crispr.find_guides(EMX1_SITE, "SpCas9")
    return next(g for g in guides if g["guide_sequence"] == "GAGTCCGAGCAGAAGAAGAA")


@pytest.fixture()
def isolated_test_ref(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIXSCOPE_REFERENCE_DIR", str(tmp_path / "refs"))
    monkeypatch.setenv("HELIXSCOPE_JOBS_DIR", str(tmp_path / "jobs"))
    genome_store.install_test_reference(fixture_path=FIXTURE)
    return genome_store.ready_record(crispr_assemblies.TEST_REFERENCE_ID)


def test_live_cas_offinder_test_reference(isolated_test_ref):
    info = crispr_casoffinder.detect_cas_offinder()
    if not info.get("available"):
        pytest.skip("Cas-OFFinder executable is not installed.")
    record = isolated_test_ref
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id=crispr_assemblies.TEST_REFERENCE_ID,
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        max_mismatches=2,
        include_nag=False,
        dna_bulge=0,
        rna_bulge=0,
        timeout_s=60,
    )
    if envelope["status"] not in {
        crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE,
        crispr_offtarget.JOB_COMPLETED,
    }:
        reason = str(envelope.get("reason") or "")
        if "OpenCL" in reason:
            pytest.skip(
                "Cas-OFFinder 2.4.1 requires a working OpenCL ICD even for "
                f"device C. This machine has no registered OpenCL device. "
                f"Search status={envelope['status']}: {reason}"
            )
        pytest.fail(f"Unexpected Cas-OFFinder status {envelope['status']}: {reason}")
    assert envelope["status"] == crispr_offtarget.JOB_COMPLETED_TEST_REFERENCE
    assert envelope["genome_wide"] is False
    assert envelope.get("tool", {}).get("available") is True
    assert envelope["verified_hit_count"] is not None
    assert envelope["verified_hit_count"] >= 1


@pytest.mark.skipif(
    os.environ.get("HELIXSCOPE_LIVE_GENOME") != "1",
    reason="Set HELIXSCOPE_LIVE_GENOME=1 only when a public assembly is READY.",
)
def test_live_public_assembly_one_guide():
    info = crispr_casoffinder.detect_cas_offinder()
    if not info.get("available"):
        pytest.skip("Cas-OFFinder executable is not installed.")
    record = genome_store.ready_record("GRCh38.p14")
    if record is None:
        pytest.skip("GRCh38.p14 is not READY on this machine.")
    envelope = crispr_casoffinder.search_indexed_fasta(
        _emx1_guide(),
        assembly_id="GRCh38.p14",
        fasta_path=record["fasta_path"],
        fai_path=record["fai_path"],
        reference_sha256=record["manifest"]["file_sha256"],
        max_mismatches=0,
        include_nag=False,
        timeout_s=3600,
        max_hits=50,
    )
    assert envelope["status"] in {
        crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE,
        crispr_offtarget.JOB_RESOURCE_LIMIT,
        crispr_offtarget.JOB_TIMEOUT,
    }
    if envelope["status"] == crispr_offtarget.JOB_COMPLETED_FULL_REFERENCE:
        assert envelope["genome_wide"] is True
    else:
        assert envelope["genome_wide"] is False
        assert envelope["guide_specificity"]["mit_sguide"]["score"] is None
