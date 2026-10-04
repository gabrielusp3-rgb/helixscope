"""Local NCBI BLAST+ live validation against a tiny declared fixture.

Remote NCBI BLAST Common URL API is a different backend and is not used here.
This module skips when BLAST+ executables are not installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from modules import blast_search

FIXTURES = Path(__file__).parent / "fixtures"
NUCL_FA = (FIXTURES / "blast_tiny_nucl.fa").read_text(encoding="utf-8")
PROT_FA = (FIXTURES / "blast_tiny_prot.fa").read_text(encoding="utf-8")
NUCL_HIT = "ATGCTAGTCGGATCCTGAATGCGTACGACTAG"
PROT_HIT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
BLASTX_DNA = (
    "ATGAAAACAGCATACATAGCAAAACAACGAACAATATCATTCGTAAAATCACACTTCTCA"
    "CGACAATTAGAAGAACGATTGGGATTAATAGAAGTACAA"
)


def _blastplus_present() -> bool:
    return bool(blast_search.local_blast_availability().get("executables_detected"))


def _makeblastdb_present() -> bool:
    return bool(blast_search.detect_makeblastdb().get("available"))


def _require_blastplus() -> None:
    if not _blastplus_present() or not _makeblastdb_present():
        pytest.skip("NCBI BLAST+ official Windows package not installed in this environment")


def test_live_makeblastdb_blastn_hit_and_nohit(tmp_path, monkeypatch) -> None:
    _require_blastplus()
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    prefix = str(tmp_path / "tiny_nucl")
    built = blast_search.build_declared_blast_database(
        fasta_text=NUCL_FA,
        dbtype="nucl",
        prefix=prefix,
        title="helixscope_tiny_nucl",
    )
    assert built["available"] is True
    assert built["not_nt"] is True
    monkeypatch.setenv("HELIXSCOPE_BLAST_DB", prefix)
    hit = blast_search.run_local_blast(
        program="blastn",
        query=NUCL_HIT,
        molecule="DNA",
    )
    assert hit["backend"] == "NCBI BLAST+ local"
    assert hit["not_remote_ncbi"] is True
    assert hit["status"] == "READY"
    assert hit["n_hits"] >= 1
    first = hit["hits"][0]
    assert "helix_nucl_hit" in str(first.get("accession") or first.get("hit_id") or "")
    assert float(first.get("identity_pct") or 0) == pytest.approx(100.0)
    assert int(first.get("alignment_length") or 0) == 32
    miss = blast_search.run_local_blast(
        program="blastn",
        query="T" * 32,
        molecule="DNA",
    )
    assert miss["backend"] == "NCBI BLAST+ local"
    assert miss["status"] == "NO_HITS"
    assert miss["n_hits"] == 0


def test_live_makeblastdb_blastp_hit(tmp_path, monkeypatch) -> None:
    _require_blastplus()
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    prefix = str(tmp_path / "tiny_prot")
    built = blast_search.build_declared_blast_database(
        fasta_text=PROT_FA,
        dbtype="prot",
        prefix=prefix,
        title="helixscope_tiny_prot",
    )
    assert built["available"] is True
    assert built["not_nr"] is True
    monkeypatch.setenv("HELIXSCOPE_BLAST_DB", prefix)
    hit = blast_search.run_local_blast(
        program="blastp",
        query=PROT_HIT,
        molecule="PROTEIN",
    )
    assert hit["backend"] == "NCBI BLAST+ local"
    assert hit["status"] == "READY"
    assert hit["n_hits"] >= 1
    first = hit["hits"][0]
    assert "helix_prot_hit" in str(first.get("accession") or first.get("hit_id") or "")
    programs = blast_search.local_blast_availability().get("programs") or {}
    for name in ("blastx", "tblastn", "tblastx"):
        row = programs.get(name) or {}
        if row.get("available"):
            assert row.get("version")


def test_live_blastx_against_tiny_prot_if_present(tmp_path, monkeypatch) -> None:
    _require_blastplus()
    programs = blast_search.local_blast_availability().get("programs") or {}
    if not (programs.get("blastx") or {}).get("available"):
        pytest.skip("blastx executable not present in this BLAST+ package")
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    prefix = str(tmp_path / "tiny_prot_x")
    blast_search.build_declared_blast_database(
        fasta_text=PROT_FA,
        dbtype="prot",
        prefix=prefix,
        title="helixscope_tiny_prot_x",
    )
    monkeypatch.setenv("HELIXSCOPE_BLAST_DB", prefix)
    hit = blast_search.run_local_blast(
        program="blastx",
        query=BLASTX_DNA,
        molecule="DNA",
    )
    assert hit["backend"] == "NCBI BLAST+ local"
    assert hit["not_remote_ncbi"] is True
    assert hit["status"] in {"READY", "NO_HITS"}


def test_live_blast_malformed_query_is_invalid_input() -> None:
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.run_local_blast(
            program="blastn",
            query="!!!!",
            molecule="DNA",
        )
    assert exc.value.category == "INVALID_INPUT"


def test_local_blast_timeout_is_timeout(monkeypatch, tmp_path) -> None:
    import subprocess

    _require_blastplus()
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    prefix = str(tmp_path / "tiny_nucl_to")
    blast_search.build_declared_blast_database(
        fasta_text=NUCL_FA,
        dbtype="nucl",
        prefix=prefix,
        title="helixscope_tiny_timeout",
    )
    monkeypatch.setenv("HELIXSCOPE_BLAST_DB", prefix)
    avail = blast_search.local_blast_availability()

    def boom(*_args, **_kwargs):
        raise subprocess.TimeoutExpired(cmd="blastn", timeout=1)

    monkeypatch.setattr(blast_search, "local_blast_availability", lambda: avail)
    monkeypatch.setattr(blast_search.subprocess, "run", boom)
    with pytest.raises(blast_search.BlastError) as exc:
        blast_search.run_local_blast(
            program="blastn",
            query=NUCL_HIT,
            molecule="DNA",
        )
    assert exc.value.category == "TIMEOUT"
    assert "NO_HITS" not in str(exc.value)
