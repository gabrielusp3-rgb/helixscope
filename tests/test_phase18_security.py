"""Auditoria Fase 18: SSRF Alignment API, caps, resource_admission, sem segredos."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from modules import rcsb_alignment, resource_admission, structure_alignment

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_alignment_api_is_in_phase18_ssrf_matrix():
    assert rcsb_alignment.url_is_allowed(
        "https://alignment.rcsb.org/api/v1/structures/submit"
    )
    assert not rcsb_alignment.url_is_allowed("https://alignment.rcsb.org/api/v1/structures/submit".replace("https", "http"))
    assert not rcsb_alignment.url_is_allowed("https://files.rcsb.org/download/1CRN.cif")


def test_resource_admission_grch38_still_resource_limit_on_this_class_of_machine():
    memory = resource_admission.system_memory()
    decision = resource_admission.admit_reference_search(
        assembly_id="GRCh38.p14",
        bases=3_100_000_000,
        engine_available=True,
    )
    assert decision["decision"] in {
        resource_admission.DECISION_RESOURCE_LIMIT,
        resource_admission.DECISION_PROCEED,
        resource_admission.DECISION_UNKNOWN,
    }
    available = int(memory.get("available_bytes") or 0)
    if memory.get("measured") and available < 6 * 1024**3:
        assert decision["decision"] == resource_admission.DECISION_RESOURCE_LIMIT
        assert decision.get("genome_wide") is not True


def test_identity_hash_includes_method_and_pair():
    payload = json.loads(
        (Path(__file__).parent / "fixtures" / "rcsb_alignment_complete.json").read_text(
            encoding="utf-8"
        )
    )
    a = structure_alignment.parse_alignment_payload(payload)
    payload["meta"]["alignment_method"] = "tm-align"
    b = structure_alignment.parse_alignment_payload(payload, method="tm-align")
    assert a["identity_hash"] != b["identity_hash"]
    assert a["method"] != b["method"]


def test_repo_has_no_hardcoded_ncbi_or_cloud_secrets():
    forbidden = re.compile(
        r"(?i)(AKIA[0-9A-Z]{16}|-----BEGIN (RSA |OPENSSH )?PRIVATE KEY-----|"
        r"sk_live_[A-Za-z0-9]+|xox[baprs]-[A-Za-z0-9-]+)"
    )
    offenders = []
    for path in REPO_ROOT.rglob("*"):
        if path.suffix.lower() not in {".py", ".md", ".json", ".txt", ".toml", ".yml", ".yaml"}:
            continue
        if any(
            part in {
                ".git",
                "node_modules",
                "__pycache__",
                "data",
                "venv",
                ".venv",
                "env",
                "ENV",
                "tools",
                "site-packages",
                ".pytest_cache",
                "htmlcov",
            }
            for part in path.parts
        ):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if forbidden.search(text):
            offenders.append(str(path.relative_to(REPO_ROOT)))
    assert offenders == []


def test_malformed_alignment_json_is_parsing_error():
    with pytest.raises(structure_alignment.StructureAlignmentError) as exc:
        structure_alignment.parse_alignment_payload({"results": "not-a-list"})
    assert exc.value.category == "PARSING_ERROR"
