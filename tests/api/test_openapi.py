"""OpenAPI contract export and stability."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from fastapi.testclient import TestClient

from helixscope_api.contract import dump_openapi, openapi_sha256
from helixscope_api.main import create_app

ROOT = Path(__file__).resolve().parents[2]
OPENAPI_PATH = ROOT / "services" / "api" / "openapi.json"

REQUIRED_PATHS = {
    "/health/live",
    "/health/ready",
    "/api/v1",
    "/api/v1/system/identity",
    "/api/v1/system/capabilities",
    "/api/v1/references",
    "/api/v1/references/admission",
    "/api/v1/references/install",
    "/api/v1/dna/analyze",
    "/api/v1/dna/helix-3d",
    "/api/v1/rna/analyze",
    "/api/v1/rna/helix-3d",
    "/api/v1/rna/translate-coding",
    "/api/v1/protein/analyze",
    "/api/v1/protein/uniprot",
    "/api/v1/protein/structure/capabilities",
    "/api/v1/protein/structure/search",
    "/api/v1/protein/structure/load",
    "/api/v1/protein/structure/analyze",
    "/api/v1/protein/structure/alphafold",
    "/api/v1/alignment/pairwise",
    "/api/v1/motif/search",
    "/api/v1/ncbi/fetch",
    "/api/v1/crispr/guides",
    "/api/v1/crispr/systems",
    "/api/v1/variants/identify",
    "/api/v1/evidence/pack",
    "/api/v1/structures/fixture",
    "/api/v1/structures/coordinate-file",
    "/api/v1/structures/experimental-complex",
    "/api/v1/compare/parse",
    "/api/v1/compare/modes",
    "/api/v1/compare/methods",
    "/api/v1/compare/variants",
    "/api/v1/compare/proteins",
    "/api/v1/compare/guides",
    "/api/v1/compare/evolution",
    "/api/v1/compare/evidence",
    "/api/v1/msa/prealigned",
    "/api/v1/phylogeny/infer",
    "/api/v1/phylogeny/taxonomy",
    "/api/v1/phylogeny/jobs",
    "/api/v1/blast/availability",
    "/api/v1/jobs/{job_id}",
}

REQUIRED_OPERATION_IDS = {
    "healthLive",
    "healthReady",
    "dnaAnalyze",
    "rnaAnalyze",
    "rnaTranslateCoding",
    "proteinAnalyze",
    "proteinUniprot",
    "proteinStructureCapabilities",
    "proteinStructureSearch",
    "proteinStructureLoad",
    "proteinStructureAnalyze",
    "proteinAlphafoldRetrieve",
    "alignmentPairwise",
    "motifSearch",
    "crisprGuides",
    "crisprListSystems",
    "variantIdentify",
    "evidencePack",
    "compareParseAlignment",
    "compareListModes",
    "compareVariants",
    "compareProteins",
    "compareGuides",
    "compareEvolution",
    "compareEvidence",
    "phylogenyInfer",
    "phylogenyAttachTaxonomy",
    "jobGet",
    "systemCapabilities",
    "systemIdentity",
    "listReferences",
    "referencesAdmitDownload",
    "referencesInstallCatalog",
    "dnaHelix3d",
    "rnaHelix3d",
    "structureLoadExperimentalComplex",
    "structureCoordinateFile",
}


def test_openapi_contract_and_export(client: TestClient) -> None:
    spec = client.app.openapi()
    dumped = json.dumps(spec, indent=2, sort_keys=True)
    assert "C:\\\\Users" not in dumped
    assert "/latest" not in spec["paths"]
    assert "streamlit" not in dumped.lower()
    paths = set(spec["paths"])
    missing = REQUIRED_PATHS - paths
    assert missing == set(), missing
    op_ids: set[str] = set()
    for item in spec["paths"].values():
        for op in item.values():
            if isinstance(op, dict) and op.get("operationId"):
                op_ids.add(str(op["operationId"]))
                assert op.get("summary")
                assert op.get("tags")
    missing_ops = REQUIRED_OPERATION_IDS - op_ids
    assert missing_ops == set(), missing_ops
    OPENAPI_PATH.write_text(dumped + "\n", encoding="utf-8", newline="\n")
    live = json.loads(OPENAPI_PATH.read_text(encoding="utf-8"))
    assert live["info"]["version"] == "v1"
    assert live["info"]["title"] == "HelixScope API"
    assert openapi_sha256(client.app) == hashlib.sha256(OPENAPI_PATH.read_bytes()).hexdigest()
    assert dump_openapi(client.app) == dumped + "\n"


def test_openapi_matches_create_app() -> None:
    spec = create_app().openapi()
    assert spec["paths"]["/api/v1/dna/analyze"]["post"]["operationId"] == "dnaAnalyze"
