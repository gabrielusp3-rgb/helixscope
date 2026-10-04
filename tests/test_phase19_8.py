"""Phase 19.8 contracts: sidebar chrome, phylogeny honesty, DSSP remote, cache."""

from __future__ import annotations

from pathlib import Path
from urllib.request import Request

import pytest

from modules import (
    dna_analysis,
    engine_validation,
    explain,
    msa,
    phylogeny,
    protein_structure,
    tool_registry,
    variant_core,
)
from ui.hand_control import HAND_CONTROL_JS

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


class _FakeHandle:
    def __init__(self, body: str, url: str) -> None:
        self._body = body.encode("utf-8")
        self._url = url

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def test_sidebar_outer_stroke_is_removed() -> None:
    styles = (ROOT / "ui" / "styles.py").read_text(encoding="utf-8")
    block = styles.split("section[data-testid=\"stSidebar\"] > div:first-child {")[1].split(
        "section[data-testid=\"stSidebar\"] > div:first-child::before"
    )[0]
    assert "border: none !important;" in block
    assert "border-right: none !important;" in block
    assert "border-left: none !important;" in block
    assert "inset 0 0 0 1px" not in block
    assert "margin: 0 6px 0 0 !important;" in block
    assert "border-radius: 0 var(--hs-radius-lg) var(--hs-radius-lg) 0 !important;" in block
    assert "backdrop-filter: blur(16px)" in block
    assert "var(--hs-glass-thick)" in block
    assert "inset 0 1px 0 var(--hs-glass-highlight)" in block
    assert "section[data-testid=\"stSidebar\"]::before" not in styles
    expand_css = styles.split("[data-testid=\"stExpandSidebarButton\"]")[1][:500]
    assert "display: none !important;" in expand_css
    shell = (ROOT / "ui" / "shell.py").read_text(encoding="utf-8")
    assert "[data-testid=\"stExpandSidebarButton\"]" in shell
    assert "[data-testid=\"stSidebarCollapseButton\"]" in shell
    collapse = shell.split("[data-testid=\"stExpandSidebarButton\"]")[1][:400]
    assert "display: none !important;" in collapse
    assert "display: flex !important;" not in collapse


def test_phylogeny_empty_copy_does_not_imply_missing_engine() -> None:
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    assert "No completed MSA is loaded." in text
    assert "Phylogeny engines are not missing" in text
    assert "Build phylogenetic tree is not shown because there is no completed MSA." not in text


def test_local_dssp_on_1crn_remains_unavailable_without_binary() -> None:
    info = protein_structure.dssp_availability()
    if info.get("available"):
        return
    assigned = protein_structure.assign_secondary_structure_dssp(
        (FIXTURES / "1CRN.cif").read_text(encoding="utf-8")
    )
    assert assigned["status"] == "UNAVAILABLE"
    assert assigned.get("engine_location") != "remote"


def test_dssp_remote_parses_injected_pdb_redo_output() -> None:
    legacy = (FIXTURES / "dssp_legacy_min.txt").read_text(encoding="utf-8")

    def fake(request: Request, timeout=None):
        assert str(request.full_url) == protein_structure.PDB_REDO_DSSP_DO
        assert request.get_method() == "POST"
        body = request.data or b""
        header = str(request.headers.get("Content-type") or request.headers.get("Content-Type") or "")
        assert "multipart/form-data" in header.lower() or b"name=\"data\"" in body
        assert b"dssp" in body
        return _FakeHandle(legacy, protein_structure.PDB_REDO_DSSP_DO)

    result = protein_structure.assign_secondary_structure_dssp_remote(
        (FIXTURES / "1CRN.cif").read_text(encoding="utf-8"),
        urlopen_fn=fake,
        coordinate_kind="experimental",
    )
    assert result["status"] == "EXPERIMENTAL"
    assert result["engine_location"] == "remote"
    assert result["tool"] == protein_structure.DSSP_REMOTE_TOOL
    assert result["not_sequence_prediction"] is True
    assert int(result.get("n_residues") or 0) >= 1
    predicted = protein_structure.assign_secondary_structure_dssp_remote(
        (FIXTURES / "1CRN.cif").read_text(encoding="utf-8")[:200],
        urlopen_fn=fake,
        coordinate_kind="predicted",
    )
    assert predicted["status"] in {"COMPUTED", "ERROR"}
    if predicted["status"] == "COMPUTED":
        assert predicted["coordinate_kind"] == "predicted"


def test_dssp_remote_empty_input_is_error_not_fake_assignment() -> None:
    result = protein_structure.assign_secondary_structure_dssp_remote("")
    assert result["status"] == "ERROR"
    assert int(result.get("n_residues") or 0) == 0


def test_structure_file_cache_roundtrip_and_tamper(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv(protein_structure.STRUCTURE_CACHE_ENV, str(tmp_path))
    url = "https://files.rcsb.org/download/1CRN.cif"
    text = "data_1CRN\n# cached authentic body\n"
    written = protein_structure.write_cached_structure_file(url, text)
    assert written["sha256"]
    hit = protein_structure.read_cached_structure_file(url)
    assert hit is not None
    assert hit[0] == text
    assert hit[1]["cache_hit"] is True
    identity = protein_structure.structure_file_cache_identity(url)
    body = tmp_path / identity[:2] / f"{identity}.cif"
    body.write_bytes(b"tampered")
    assert protein_structure.read_cached_structure_file(url) is None
    assert protein_structure.read_cached_structure_file("https://evil.example/1CRN.cif") is None


def test_tool_registry_keeps_local_and_remote_dssp_distinct() -> None:
    snap = tool_registry.collect_tool_snapshot(refresh=True)
    tools = snap["tools"]
    assert "DSSP" in tools
    assert "DSSP PDB-REDO API" in tools
    assert tools["DSSP PDB-REDO API"]["engine_location"] == "remote"
    assert tools["DSSP"].get("engine_location") != "remote"
    remote_live = tools["DSSP PDB-REDO API"].get("live") or {}
    if remote_live:
        assert remote_live.get("status") != engine_validation.STATUS_LIVE_VALIDATED
        assert remote_live.get("status") in {
            engine_validation.STATUS_REMOTE_VALIDATED,
            engine_validation.STATUS_BROKEN,
            engine_validation.STATUS_DETECTED,
        } or remote_live.get("engine_location") == "remote"


def test_cached_snapshot_overlays_live_record(monkeypatch) -> None:
    tool_registry.collect_tool_snapshot(refresh=True)
    real = engine_validation.live_record

    def fake_live(name: str):
        if name == "IQ-TREE":
            return {
                "tool": "IQ-TREE",
                "ok": True,
                "status": engine_validation.STATUS_LIVE_VALIDATED,
                "version": "overlay-test",
            }
        return real(name)

    monkeypatch.setattr(engine_validation, "live_record", fake_live)
    snap = tool_registry.collect_tool_snapshot(refresh=False)
    assert snap["tools"]["IQ-TREE"]["live"]["version"] == "overlay-test"
    assert snap["tools"]["IQ-TREE"]["live"]["status"] == engine_validation.STATUS_LIVE_VALIDATED


def test_nj_and_upgma_consume_the_same_msa_without_invented_support() -> None:
    aligned = msa.import_prealigned_fasta(
        ">seq_a\nACGTACGTACGTACGTACGT\n"
        ">seq_b\nACGTACGTACGTACGTTTTT\n"
        ">seq_c\nTTTTACGTACGTACGTACGT\n"
        ">seq_d\nACGTACGTAAAAACGTACGT\n"
    )
    nj = phylogeny.infer_phylogeny(
        aligned,
        method=phylogeny.METHOD_NJ,
        rooting=phylogeny.ROOTING_UNROOTED,
    )
    upgma = phylogeny.infer_phylogeny(
        aligned,
        method=phylogeny.METHOD_UPGMA,
    )
    assert nj["status"] == "COMPUTED"
    assert upgma["status"] == "COMPUTED"
    assert nj["n_leaves"] == 4
    assert upgma["n_leaves"] == 4
    assert nj["newick"].count("seq_") == 4
    assert upgma["newick"].count("seq_") == 4
    assert nj["alignment_hash"] == upgma["alignment_hash"]
    assert nj["alignment_hash"]
    for tree in (nj, upgma):
        support = tree.get("support") or {}
        methods = [str(item.get("method")) for item in list(support.get("methods") or [])]
        assert phylogeny.SUPPORT_UFBOOT not in methods
        assert phylogeny.SUPPORT_SH_ALRT not in methods


def test_hand_control_plotly_lookup_uses_public_widget_key() -> None:
    assert "st-key-" in HAND_CONTROL_JS
    assert "plotHostKey" in HAND_CONTROL_JS
    assert "js-plotly-plot" in HAND_CONTROL_JS
    source = (ROOT / "app.py").read_text(encoding="utf-8")
    assert 'st.container(key=f"{key}_host")' in source


def test_explain_does_not_equate_conservation_with_pathogenicity() -> None:
    rec = explain.explain(
        "evolution",
        {"status": "COMPUTED", "conservation_label": "high"},
    )
    blob = f"{rec.plain_meaning} {rec.limitations} {rec.interpretation}".lower()
    assert "not prove" in blob or "does not prove" in blob
    assert "pathogenic" in blob
    assert "conservation is not function" in blob


def test_reverse_complement_involution_ascii_dna() -> None:
    seq = "ATGCATGC"
    twice = dna_analysis.reverse_complement(dna_analysis.reverse_complement(seq))
    assert twice == seq


def test_adversarial_newick_is_parsing_error() -> None:
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.parse_newick("); DROP TABLE trees; --")
    assert exc.value.category == "PARSING_ERROR"


def test_parse_ply_ascii_rejects_empty_and_accepts_triangle() -> None:
    with pytest.raises(protein_structure.StructureError) as exc:
        protein_structure.parse_ply_ascii("")
    assert exc.value.category == "PARSING_ERROR"
    ply = (
        "ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\n"
        "property float y\nproperty float z\nelement face 1\n"
        "property list uchar int vertex_indices\nend_header\n"
        "0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n"
    )
    mesh = protein_structure.parse_ply_ascii(ply)
    assert mesh["n_vertices"] == 3
    assert mesh["n_faces"] == 1
    assert mesh["faces"] == [[0, 1, 2]]


def test_variant_parser_rejects_markup_and_path_traversal() -> None:
    with pytest.raises(variant_core.VariantInputError):
        variant_core.parse_variant_input("<script>alert(1)</script>")
    with pytest.raises(variant_core.VariantInputError):
        variant_core.parse_variant_input("../../etc/passwd")
    with pytest.raises(variant_core.VariantInputError):
        variant_core.parse_variant_input("NM_000000.0:c.1A>T")


def test_azimuth_and_deephf_remain_unavailable_not_zero() -> None:
    from modules import crispr_ontarget

    rs2 = crispr_ontarget.ruleset2_availability()
    deephf = crispr_ontarget.deephf_availability()
    assert rs2["available"] is False
    assert deephf["available"] is False
    scored = crispr_ontarget.score_doench_ruleset2("ACGT" * 8)
    assert scored.get("score") is None
    assert scored.get("status") == "UNAVAILABLE"
