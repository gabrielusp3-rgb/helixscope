"""Testes de admissao de recursos. Sem rede, sem tarefas pesadas reais."""

from __future__ import annotations

import pytest

from modules import resource_admission

GIB = 1024**3


def _memory(total_gib: float, available_gib: float, page_gib: float = 0.0) -> dict:
    return {
        "measured": True,
        "source": "injected",
        "total_bytes": int(total_gib * GIB),
        "available_bytes": int(available_gib * GIB),
        "load_percent": int(round(100 * (total_gib - available_gib) / total_gib)),
        "page_file_total_bytes": int(page_gib * GIB),
        "page_file_available_bytes": int(page_gib * GIB),
    }


def test_system_memory_reports_a_measurement_or_says_it_cannot():
    reading = resource_admission.system_memory()
    assert set(reading) >= {"measured", "source", "total_bytes", "available_bytes"}
    if reading["measured"]:
        assert reading["total_bytes"] > 0
        assert reading["available_bytes"] >= 0
        assert reading["available_bytes"] <= reading["total_bytes"]
        assert reading["source"] != resource_admission.MEMORY_SOURCE_UNAVAILABLE
    else:
        assert reading["total_bytes"] is None
        assert reading["reason"]


def test_disk_free_measures_the_repository_volume(tmp_path):
    disk = resource_admission.disk_free(str(tmp_path))
    assert disk["measured"] is True
    assert disk["total_bytes"] > 0
    assert 0 <= disk["free_bytes"] <= disk["total_bytes"]


def test_disk_free_walks_up_to_an_existing_parent(tmp_path):
    disk = resource_admission.disk_free(str(tmp_path / "does" / "not" / "exist.fa"))
    assert disk["measured"] is True


def test_base_count_estimate_discounts_headers_and_newlines():
    bases = resource_admission.estimate_bases_from_fasta_bytes(1_000_000)
    assert bases < 1_000_000
    assert bases == int(1_000_000 * (1 - resource_admission.FASTA_HEADER_AND_NEWLINE_OVERHEAD))
    assert resource_admission.estimate_bases_from_fasta_bytes(0) == 0
    with pytest.raises(ValueError):
        resource_admission.estimate_bases_from_fasta_bytes(-1)


def test_memory_estimate_is_labelled_as_an_estimate_with_a_declared_basis():
    estimate = resource_admission.estimate_cas_offinder_memory(bases=3_100_000_000)
    assert estimate["is_estimate"] is True
    assert estimate["bases"] == 3_100_000_000
    assert estimate["estimated_bytes"] > 3_100_000_000
    assert "one byte per base" in estimate["basis"]
    assert "not a measurement" in estimate["basis"]
    assert estimate["engine"] == "Cas-OFFinder 2.4.1"
    from_bytes = resource_admission.estimate_cas_offinder_memory(fasta_bytes=3_200_000_000)
    assert from_bytes["bases"] < 3_200_000_000
    with pytest.raises(ValueError):
        resource_admission.estimate_cas_offinder_memory()


def test_low_memory_machine_gets_resource_limit_not_a_partial_score():
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        memory=_memory(total_gib=7.8, available_gib=1.2),
    )
    assert decision["decision"] == resource_admission.DECISION_RESOURCE_LIMIT
    assert decision["blocker"] == "memory"
    assert "RESOURCE_LIMIT is the correct result" in decision["reason"]
    assert "no genome-wide specificity score is produced" in decision["reason"]


def test_machine_with_room_is_admitted():
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        memory=_memory(total_gib=32.0, available_gib=24.0),
    )
    assert decision["decision"] == resource_admission.DECISION_PROCEED
    assert decision["blocker"] == ""
    assert decision["budget_bytes"] >= decision["required_bytes"]


def test_page_file_is_not_counted_unless_explicitly_allowed():
    memory = _memory(total_gib=7.8, available_gib=1.2, page_gib=16.0)
    refused = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14", bases=3_100_000_000, memory=memory
    )
    assert refused["decision"] == resource_admission.DECISION_RESOURCE_LIMIT
    assert refused["page_file_counted"] is False
    allowed = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        memory=memory,
        allow_page_file=True,
    )
    assert allowed["decision"] == resource_admission.DECISION_PROCEED
    assert allowed["page_file_counted"] is True


def test_unmeasurable_memory_is_unknown_and_never_proceed():
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        memory={"measured": False, "available_bytes": None, "total_bytes": None},
    )
    assert decision["decision"] == resource_admission.DECISION_UNKNOWN
    assert decision["decision"] != resource_admission.DECISION_PROCEED
    assert "UNKNOWN is not PROCEED" in decision["reason"]


def test_missing_engine_blocks_before_any_memory_arithmetic():
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        engine_available=False,
        engine_reason="Cas-OFFinder is NOT_INSTALLED.",
    )
    assert decision["decision"] == resource_admission.DECISION_RESOURCE_LIMIT
    assert decision["blocker"] == "engine"
    assert decision["estimate"] is None
    assert "NOT_INSTALLED" in decision["reason"]


def test_download_admission_counts_archive_and_decompressed_copies(tmp_path):
    decision = resource_admission.admit_download(
        assembly_id="GRCh38.p14",
        compressed_bytes=1024**3,
        store_path=str(tmp_path),
    )
    assert decision["decision"] in {
        resource_admission.DECISION_PROCEED,
        resource_admission.DECISION_RESOURCE_LIMIT,
    }
    assert decision["required_bytes"] > 3 * 1024**3
    assert decision["gzip_ratio_assumed"] == 3.0
    without_archive = resource_admission.admit_download(
        assembly_id="GRCh38.p14",
        compressed_bytes=1024**3,
        store_path=str(tmp_path),
        keep_archive=False,
    )
    assert without_archive["required_bytes"] < decision["required_bytes"]
    with pytest.raises(ValueError):
        resource_admission.admit_download(
            assembly_id="GRCh38.p14", compressed_bytes=0, store_path=str(tmp_path)
        )


def test_decision_records_are_traceable():
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=1000,
        memory=_memory(total_gib=32.0, available_gib=24.0),
    )
    assert decision["assembly_id"] == "GRCh38.p14"
    assert decision["checked_at_utc"].endswith("Z")
    assert decision["software_version"]
    assert decision["estimate"]["overhead_factor"] == resource_admission.CAS_OFFINDER_OVERHEAD_FACTOR
