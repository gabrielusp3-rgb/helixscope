"""Testes de estruturas proteicas: PDB/mmCIF, mapping, fontes e seguranca.

1CRN.cif e rcsb_entry_1CRN.json sao copias publicas do RCSB PDB (crambin,
X-ray, 1.5 A). alphafold_prediction_P01308.json e metadados publicos da
AlphaFold DB para INS_HUMAN (P01308). uniprot_P01308_structure_xrefs.json e
um recorte UniProtKB. Fixtures parser_* testam o parser e nao sao deposicoes.
Nenhum teste inventa coordenadas, pLDDT ou resolucao.
"""

from __future__ import annotations

import io
import json
from email.message import EmailMessage
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request

import pytest

from modules import protein_structure, provenance, scientific_checks, tool_detection

FIXTURES = Path(__file__).parent / "fixtures"
CRAMBIN = "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
EMAIL_NOT_USED = "unused"


class _FakeHandle:
    def __init__(self, body: str | bytes, url: str = "") -> None:
        self._body = body.encode("utf-8") if isinstance(body, str) else body
        self._url = url

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


def _load(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def test_identifier_classification_rejects_urls_and_accepts_real_ids():
    assert protein_structure.classify_structure_identifier("1crn")["kind"] == "pdb"
    assert protein_structure.classify_structure_identifier("P01308")["kind"] == "uniprot"
    af = protein_structure.classify_structure_identifier("AF-P01308-F1")
    assert af["kind"] == "alphafold"
    assert af["uniprot"] == "P01308"
    assert protein_structure.classify_structure_identifier("https://files.rcsb.org/download/1CRN.cif")["kind"] == "invalid"
    assert protein_structure.classify_structure_identifier("../etc/passwd")["kind"] == "invalid"


def test_url_allowlist_blocks_ssrf():
    assert protein_structure.url_is_allowed("https://files.rcsb.org/download/1CRN.cif")
    assert protein_structure.url_is_allowed("https://data.rcsb.org/rest/v1/core/entry/1CRN")
    assert protein_structure.url_is_allowed("https://alphafold.ebi.ac.uk/api/prediction/P01308")
    assert protein_structure.url_is_allowed("https://rest.uniprot.org/uniprotkb/P01308.json") is True
    assert protein_structure.url_is_allowed("https://pdb-redo.eu/dssp/do") is True
    assert protein_structure.url_is_allowed("https://pdb-redo.eu/dssp/db/1crn/legacy") is True
    assert protein_structure.url_is_allowed("https://pdb-redo.eu/dssp/about") is False
    assert protein_structure.url_is_allowed("http://pdb-redo.eu/dssp/do") is False
    assert protein_structure.url_is_allowed("http://files.rcsb.org/download/1CRN.cif") is False
    assert protein_structure.url_is_allowed("https://127.0.0.1/1CRN.cif") is False
    assert protein_structure.url_is_allowed("https://evil.example/1CRN.cif") is False
    with pytest.raises(protein_structure.StructureError) as exc:
        protein_structure.fetch_structure_file("https://evil.example/1CRN.cif")
    assert exc.value.category == "INVALID_INPUT"


def test_parse_real_1crn_mmcif_and_rcsb_metadata():
    parsed = protein_structure.parse_mmcif(_load("1CRN.cif"))
    assert parsed["entry_id"] == "1CRN"
    assert parsed["method"] == "X-RAY DIFFRACTION"
    assert parsed["resolution_angstrom"] == pytest.approx(1.5)
    assert parsed["n_atoms"] == 327
    chain = parsed["chains"][0]
    assert chain["sequence"] == CRAMBIN
    assert chain["chain_id"] == "A"
    assert sum(1 for item in chain["residues"] if item["has_coordinates"]) == 46
    meta = protein_structure.parse_rcsb_entry_metadata(json.loads(_load("rcsb_entry_1CRN.json")))
    assert meta["kind"] == "experimental"
    assert meta["source"] == "RCSB PDB"
    assert meta["resolution_angstrom"] == pytest.approx(1.5)
    assert meta["method"] == "X-RAY DIFFRACTION"


def test_1crn_mapping_is_not_auth_equals_index():
    parsed = protein_structure.parse_mmcif(_load("1CRN.cif"))
    mapped = protein_structure.map_query_to_chain(CRAMBIN, parsed["chains"][0])
    assert mapped["sequences_identical"] is True
    assert mapped["coverage_query_with_coordinates"] == pytest.approx(1.0)
    assert scientific_checks.protein_mapping_covers_query(len(mapped["mapping"]), 46)
    first = mapped["mapping"][0]
    assert first["query_index_0based"] == 0
    assert first["label_seq_id"] == 1
    assert first["auth_seq_id"] == 1
    truncated = protein_structure.map_query_to_chain(CRAMBIN[:20], parsed["chains"][0])
    assert truncated["sequences_identical"] is False
    assert truncated["n_query"] == 20
    mutated = "A" + CRAMBIN[1:]
    mismatch = protein_structure.map_query_to_chain(mutated, parsed["chains"][0])
    assert mismatch["mapping"][0]["match"] is False
    assert mismatch["identity_ungapped_pct"] < 100.0
    assert mapped["mapping_status"] == "EXACT"
    assert scientific_checks.mapping_status_is_declared(mapped["mapping_status"])
    assert mismatch["mapping_status"] == "ALIGNED"
    assert "not exact" in mismatch["mapping_status_note"].lower()


def test_nmr_models_altloc_insertion_and_missing_residue():
    parsed = protein_structure.parse_mmcif(_load("parser_nmr_altloc.cif"))
    assert parsed["models"] == [1, 2]
    assert parsed["method"] == "NMR"
    model_one = next(item for item in parsed["chains"] if item["model"] == 1)
    by_label = {item["label_seq_id"]: item for item in model_one["residues"]}
    assert by_label[2]["has_coordinates"] is False
    assert by_label[1]["alt_ids"] == ["A", "B"]
    assert by_label[3]["insertion_code"] == "A"
    assert by_label[1]["auth_seq_id"] == 10
    assert by_label[1]["auth_seq_id"] != by_label[1]["label_seq_id"]
    mapped = protein_structure.map_query_to_chain("AKA", model_one)
    assert mapped["mapping"][1]["has_coordinates"] is False
    assert mapped["mapping"][1]["coordinate_status"] == "unavailable"
    assert mapped["mapping"][1]["sequence_residue_exists"] is True
    assert mapped["n_without_coordinates"] == 1
    assert mapped["mapping_status"] == "EXACT"
    assert mapped["mapping"][1]["mapping_row_status"] == "EXACT"
    assert mapped["mapping"][1]["sequence_residue_exists"] is True


def test_legacy_pdb_parser_and_malformed_coordinates():
    parsed = protein_structure.parse_pdb(_load("parser_legacy.pdb"))
    assert parsed["n_atoms"] == 2
    assert parsed["chains"][0]["sequence"] == "AG"
    with pytest.raises(protein_structure.StructureError) as exc:
        protein_structure.parse_pdb("ATOM      1  CA  ALA A   1       nan   2.000   3.000")
    assert exc.value.category == "PARSING_ERROR"
    with pytest.raises(protein_structure.StructureError) as empty:
        protein_structure.parse_mmcif("   ")
    assert empty.value.category == "PARSING_ERROR"


def test_pdb_auth_numbering_is_not_treated_as_polymer_index():
    parsed = protein_structure.parse_pdb(_load("parser_auth_offset.pdb"))
    chain = parsed["chains"][0]
    assert chain["sequence"] == "AK"
    assert len(chain["residues"]) == 2
    labels = {item["label_seq_id"] for item in chain["residues"]}
    assert 1 not in labels
    assert 10 in labels
    mapped = protein_structure.map_query_to_chain("AK", chain)
    assert mapped["mapping"][0]["query_index_0based"] == 0
    assert mapped["mapping"][0]["auth_seq_id"] == 10
    assert mapped["mapping"][0]["label_seq_id"] == 10
    assert mapped["mapping"][0]["has_coordinates"] is True
    assert mapped["mapping"][1]["auth_seq_id"] == 11
    assert mapped["coverage_query_with_coordinates"] == pytest.approx(1.0)
    assert mapped["mapping_status"] == "PARTIAL"
    assert mapped["mapping_status"] != "EXACT"
    assert "not exact" in mapped["mapping_status_note"].lower()


def test_alphafold_metadata_is_predicted_not_experimental():
    payload = json.loads(_load("alphafold_prediction_P01308.json"))
    parsed = protein_structure.parse_alphafold_prediction(payload)
    assert parsed["kind"] == "predicted"
    assert parsed["source"] == "AlphaFold DB"
    assert parsed["uniprot"] == "P01308"
    assert parsed["global_metric_value"] == pytest.approx(52.91)
    assert "probability" not in parsed["global_metric_name"].lower() or "not a probability" in parsed["global_metric_name"].lower()
    text = protein_structure.source_disclaimer("predicted", "AlphaFold DB")
    assert "Predicted" in text
    assert "AlphaFold DB" in text
    assert "experimental measurement" in text.lower() or "not an experimental" in text.lower()


def test_uniprot_xrefs_are_truncated_not_invented():
    payload = json.loads(_load("uniprot_P01308_structure_xrefs.json"))
    parsed = protein_structure.parse_uniprot_structure_xrefs(payload)
    assert parsed["accession"] == "P01308"
    assert parsed["n_pdb_reported"] >= len(parsed["pdb"])
    assert parsed["truncated"] is True
    assert parsed["pdb"][0]["kind"] == "experimental"
    assert parsed["pdb"][0]["match_basis"].startswith("UniProtKB")
    nmr = next(item for item in parsed["pdb"] if item["structure_id"] == "1A7F")
    assert nmr["resolution_angstrom"] is None
    assert nmr["method"] == "NMR"


def test_load_1crn_without_network_and_cache():
    meta = protein_structure.parse_rcsb_entry_metadata(json.loads(_load("rcsb_entry_1CRN.json")))
    cache: dict = {}
    result = protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        cache=cache,
        urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
    )
    assert result["status"] == "RETRIEVED"
    assert result["kind"] == "experimental"
    assert result["source"] == "RCSB PDB"
    assert result["structure_id"] == "1CRN"
    assert result["sequence_hash"] == provenance.sequence_digest(CRAMBIN)
    assert result["cache_status"] == "live"
    assert result["resolution_angstrom"] == pytest.approx(1.5)
    assert result["coverage_query_with_coordinates"] == pytest.approx(1.0)
    assert "experimental" in result["disclaimer"].lower()
    pair = protein_structure.query_position_to_residue(result, 0)
    assert pair["auth_seq_id"] == 1
    assert protein_structure.map_msa_column_to_structure_residue([None, 0, 1], 0, result) is None
    assert protein_structure.map_msa_column_to_structure_residue([None, 0, 1], 2, result)["query_index_0based"] == 1
    span = protein_structure.residues_for_sequence_span(result, 0, 2)
    assert len(span) == 2
    again = protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        cache=cache,
        urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
    )
    assert again["cache_status"] == "cached"
    assert again["retrieved_at"] == result["retrieved_at"]
    other = protein_structure.load_protein_structure(
        CRAMBIN[:20],
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        cache=cache,
        urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
    )
    assert other["cache_status"] == "live"
    assert other["sequence_hash"] != result["sequence_hash"]
    bundle = protein_structure.export_result_bundle(result)
    assert "raw_structure_text" not in bundle
    assert bundle["resolution_angstrom"] == pytest.approx(1.5)
    csv_text = protein_structure.export_mapping_csv(result)
    assert "1CRN" in csv_text
    assert "N/A" not in csv_text.splitlines()[1].split(",")[0]
    view = protein_structure.structure_3d_input(result)
    assert view["renderer"] == "data-only"
    assert view["kind"] == "experimental"
    assert view["mapping_status"] == "EXACT"
    reverse = protein_structure.structure_residue_to_query(
        result, chain_id="A", label_seq_id=1, auth_seq_id=1
    )
    assert reverse["query_index_0based"] == 0
    assert protein_structure.structure_residue_to_query(
        result, chain_id="A", auth_seq_id=999
    ) is None
    assert view["sequence_hash"] == result["sequence_hash"]
    assert view["structure_hash"] == result["content_hash"]
    assert len(view["coordinates"]) == 327
    dumped = str(view).lower()
    assert "http://" not in dumped
    assert "https://" not in dumped
    original_x = result["atoms"][0]["x"]
    view["coordinates"][0]["x"] = 999.0
    view["atoms"][0]["x"] = 999.0
    assert result["atoms"][0]["x"] == original_x
    assert result["sequence_hash"] == view["sequence_hash"]


def test_search_sources_mocked_and_no_hits_is_not_a_structure():
    entry = _load("rcsb_entry_1CRN.json")
    search = _load("rcsb_search_1CRN.json")

    def fake(request: Request, timeout=None):
        url = str(request.full_url)
        if "core/entry/1CRN" in url:
            return _FakeHandle(entry, url)
        if "rcsbsearch" in url:
            return _FakeHandle(search, url)
        raise AssertionError(url)

    found = protein_structure.search_structure_sources(
        CRAMBIN,
        pdb_id="1CRN",
        search_pdb_by_sequence=True,
        urlopen_fn=fake,
    )
    assert found["status"] == "RETRIEVED"
    assert found["hits"]
    assert found["hits"][0]["kind"] == "experimental"
    empty = protein_structure.search_structure_sources(
        CRAMBIN,
        search_pdb_by_sequence=True,
        urlopen_fn=lambda request, timeout=None: _FakeHandle("", str(request.full_url)),
    )
    assert empty["status"] == "UNAVAILABLE"
    assert empty["category"] == "NO_STRUCTURE"
    assert empty["hits"] == []


def test_http_errors_are_classified(monkeypatch):
    def boom_404(request, timeout=None):
        raise HTTPError(
            str(request.full_url),
            404,
            "missing",
            hdrs=EmailMessage(),
            fp=io.BytesIO(b""),
        )

    with pytest.raises(protein_structure.StructureError) as missing:
        protein_structure.fetch_rcsb_entry("1CRN", urlopen_fn=boom_404)
    assert missing.value.category == "NOT_FOUND"

    def boom_timeout(request, timeout=None):
        raise TimeoutError("timed out")

    with pytest.raises(protein_structure.StructureError) as timed:
        protein_structure.fetch_rcsb_entry("1CRN", urlopen_fn=boom_timeout)
    assert timed.value.category == "TIMEOUT"

    def boom_big(request, timeout=None):
        return _FakeHandle("A" * (protein_structure.MAX_METADATA_BYTES + 1), str(request.full_url))

    with pytest.raises(protein_structure.StructureError) as limited:
        protein_structure.fetch_rcsb_entry("1CRN", urlopen_fn=boom_big)
    assert limited.value.category == "RESOURCE_LIMIT"

    def boom_net(request, timeout=None):
        raise URLError("connection refused")

    with pytest.raises(protein_structure.StructureError) as net:
        protein_structure.fetch_rcsb_entry("1CRN", urlopen_fn=boom_net)
    assert net.value.category == "NETWORK_ERROR"
    assert net.value.category != "NOT_FOUND"
    assert net.value.category != "NO_STRUCTURE"

    def boom_json(request, timeout=None):
        return _FakeHandle("{not json", str(request.full_url))

    with pytest.raises(protein_structure.StructureError) as parse:
        protein_structure.fetch_rcsb_entry("1CRN", urlopen_fn=boom_json)
    assert parse.value.category == "PARSING_ERROR"


def test_bundled_1crn_hit_matches_crambin_without_network() -> None:
    hit = protein_structure.bundled_deposited_hit("1CRN", CRAMBIN, "PROTEIN")
    assert hit is not None
    assert hit["structure_id"] == "1CRN"
    assert hit["kind"] == "experimental"
    assert hit["bundled_fixture"] == "1CRN.cif"
    assert "live" not in str(hit.get("match_basis") or "").lower() or "not a live" in str(hit.get("match_basis") or "").lower()
    missing = protein_structure.bundled_deposited_hit("1CRN", "ACDEFGHIKLMNPQRSTVWY", "PROTEIN")
    assert missing is None
    dna = protein_structure.bundled_deposited_hit("1BNA", "CGCGAATTCGCG", "DNA")
    assert dna is not None
    assert dna["bundled_fixture"] == "1BNA.cif"
    text = protein_structure.read_bundled_mmcif("1CRN.cif")
    assert "_entry.id" in text
    with pytest.raises(protein_structure.StructureError) as bad:
        protein_structure.read_bundled_mmcif("../1CRN.cif")
    assert bad.value.category == "INVALID_INPUT"


def test_oversized_structure_is_resource_limit_not_empty_coords(monkeypatch):
    huge = "data_X\nloop_\n_atom_site.Cartn_x\n" + ("1.0\n" * 10) + "extra"
    monkeypatch.setattr(protein_structure, "MAX_STRUCTURE_BYTES", 20)
    with pytest.raises(protein_structure.StructureError) as exc:
        protein_structure.parse_mmcif(huge)
    assert exc.value.category == "RESOURCE_LIMIT"


def test_coverage_zero_length_is_none_not_zero_fake():
    assert protein_structure.coverage_from_mapping(0, 0) is None
    assert protein_structure.coverage_from_mapping(10, 0) == pytest.approx(0.0)
    csv_text = protein_structure.export_mapping_csv(
        {
            "residue_mapping": [],
            "coverage_query_with_coordinates": 0.0,
            "identity_ungapped_pct": float("nan"),
            "resolution_angstrom": None,
            "structure_id": "1CRN",
            "source": "RCSB PDB",
            "kind": "experimental",
            "method": "X-RAY DIFFRACTION",
            "sequence_hash": "abc",
            "status": "RETRIEVED",
        }
    )
    assert "0.0" in csv_text or ",0," in csv_text
    assert "N/A" in csv_text


def test_best_effort_mapping_is_not_called_exact():
    parsed = protein_structure.parse_mmcif(_load("parser_best_effort.cif"))
    chain = parsed["chains"][0]
    assert chain["sequence"] == "AGKA"
    mapped = protein_structure.map_query_to_chain("AGKA", chain)
    assert mapped["mapping_status"] == "BEST_EFFORT"
    assert mapped["mapping_status"] != "EXACT"
    assert "not exact" in mapped["mapping_status_note"].lower()
    assert mapped["mapping"][0]["auth_seq_id"] == 10
    assert mapped["mapping"][0]["query_index_0based"] == 0
    assert mapped["mapping"][1]["has_coordinates"] is False
    assert mapped["mapping"][2]["auth_seq_id"] == 20


def test_remote_source_status_reports_plotly_renderer():
    info = protein_structure.remote_source_status()
    assert info["rcsb_pdb"]["available"] is True
    assert info["alphafold_db"]["kind"] == "predicted"
    assert info["local_3d_renderer"]["available"] is True
    assert info["local_3d_renderer"]["engine"] == "plotly"


def test_parse_1crn_completes_within_two_seconds():
    import time

    text = _load("1CRN.cif")
    started = time.perf_counter()
    protein_structure.parse_mmcif(text)
    elapsed = time.perf_counter() - started
    assert elapsed < 2.0


def test_json_export_keeps_none_for_missing_resolution():
    sanitized = provenance.json_safe({"resolution_angstrom": float("nan"), "n_atoms": 0})
    assert sanitized["resolution_angstrom"] is None
    assert sanitized["n_atoms"] == 0


def test_environment_matrix_is_not_invented():
    report = tool_detection.environment_report(tools={"RNAfold": {"available": False, "version": "", "path": ""}})
    assert report["python_version"].startswith("3.")
    assert "tools" in report


def test_protein_structure_section_is_before_local_ss_unavailable():
    source = Path(__file__).resolve().parents[1] / "app.py"
    text = source.read_text(encoding="utf-8")
    start = text.index("def render_protein_analysis")
    body = text[start:]
    assert body.index("_render_protein_structure_section") < body.index('key="prot_secondary"')
    assert "scatter_3d" not in body
    assert "View 3D Structure" in text
    assert "Structure retrieval failed." in text


def test_protein_analyze_opens_structure_section_without_fake_coords():
    from tests.helix_apptest import open_module

    app = open_module("protein", timeout=60)
    assert not app.exception
    app.text_area(key="prot_text").set_value(CRAMBIN).run()
    assert not app.exception
    app.button(key="prot_analyze").click().run()
    assert not app.exception
    app.toggle(key="prot_structure").set_value(True).run()
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("info", "error", "markdown", "caption", "warning")
        for item in list(getattr(app, bucket, []) or [])
    )
    keys = [item.key for item in app.button]
    assert "prot_struct_search_btn" in keys
    assert "prot_struct_load" not in keys
    toggle_keys = [item.key for item in app.toggle]
    assert "prot_view_3d" not in toggle_keys
    lowered = combined.lower()
    assert "view 3d structure" in lowered or "plotly" in lowered
    assert "does not imply a structure" in lowered or "not imply a structure" in lowered


def test_view_3d_toggle_renders_injected_1crn_without_network():
    from tests.helix_apptest import open_module

    meta = protein_structure.parse_rcsb_entry_metadata(json.loads(_load("rcsb_entry_1CRN.json")))
    loaded = protein_structure.load_protein_structure(
        CRAMBIN,
        source="RCSB PDB",
        structure_id="1CRN",
        structure_text=_load("1CRN.cif"),
        metadata=meta,
        urlopen_fn=lambda *_a, **_k: (_ for _ in ()).throw(AssertionError("network")),
    )
    app = open_module("protein", timeout=60)
    app.text_area(key="prot_text").set_value(CRAMBIN).run()
    app.button(key="prot_analyze").click().run()
    app.session_state["prot_struct_result"] = loaded
    app.toggle(key="prot_structure").set_value(True).run()
    assert not app.exception
    toggle_keys = [item.key for item in app.toggle]
    assert "prot_view_3d" in toggle_keys
    app.toggle(key="prot_view_3d").set_value(True).run()
    assert not app.exception
    combined = " ".join(
        str(getattr(item, "value", item))
        for bucket in ("info", "error", "markdown", "caption", "warning")
        for item in list(getattr(app, bucket, []) or [])
    ).lower()
    assert "experimental" in combined
    assert "1crn" in combined
    assert "structure inspector" in combined
    assert "internal position" in combined
    assert "partial structure view" not in combined or "selected chain" in combined
    assert "http://" not in combined
    assert "molecular dynamics" not in combined or "not molecular dynamics" in combined
    errors = [str(item.value) for item in list(app.error) or []]
    assert errors == []
    button_keys = [item.key for item in app.button]
    assert "msa_3d_highlight" not in button_keys
    assert "motif_3d_highlight" not in button_keys

