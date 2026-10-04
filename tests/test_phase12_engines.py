"""Testes da Fase 12: ML IQ-TREE, DSSP mapping, SASA, CRISPR models, BLAST+, registry."""

from __future__ import annotations

import inspect
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from modules import (
    blast_search,
    crispr,
    crispr_ontarget,
    msa,
    phylogeny,
    protein_structure,
    tool_registry,
)

FIXTURES = Path(__file__).parent / "fixtures"


def _msa_from_aligned(
    identifiers: list[str],
    aligned: list[str],
    *,
    molecule: str = "DNA",
) -> dict:
    members = []
    raw_lines = []
    for ident, row in zip(identifiers, aligned):
        ungapped = row.replace("-", "").replace(".", "")
        members.append(
            msa.build_collection_member(
                sequence=ungapped,
                identifier=ident,
                source="user input",
                molecule=molecule,
                organism="",
                accession="",
                version="",
            )
        )
        raw_lines.append(f">{ident}")
        raw_lines.append(row)
    return msa.build_msa_result(
        raw_alignment="\n".join(raw_lines) + "\n",
        members=members,
        tool="fixture",
        tool_version="test",
        method="prealigned fixture",
        parameters={"source": "test fixture"},
    )


def test_iqtree_absent_is_unavailable(monkeypatch):
    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: None,
    )
    info = phylogeny.detect_iqtree()
    assert info["available"] is False
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(result, method=phylogeny.METHOD_IQTREE)
    assert exc.value.category == "TOOL_NOT_INSTALLED"


def test_iqtree_ufboot_below_minimum_is_invalid():
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(
            result,
            method=phylogeny.METHOD_IQTREE,
            ufboot_replicates=100,
        )
    assert exc.value.category == "INVALID_INPUT"


def test_iqtree_modelfinder_and_support_from_mocked_binary(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    validated = phylogeny.validate_msa_for_phylogeny(result)
    names = [item["tree_id"] for item in validated["leaf_records"]]
    newick = f"({names[0]}:0.1,({names[1]}:0.2,{names[2]}:0.3)95.4/98:0.05);"

    def fake_run(argv, **kwargs):
        assert kwargs.get("shell") is False
        assert "-m" in argv and argv[argv.index("-m") + 1] == "MFP"
        assert "-B" in argv and argv[argv.index("-B") + 1] == "1000"
        assert "-alrt" in argv and argv[argv.index("-alrt") + 1] == "1000"
        prefix = argv[argv.index("--prefix") + 1]
        Path(prefix + ".treefile").write_text(newick, encoding="utf-8")
        Path(prefix + ".iqtree").write_text(
            "Best-fit model according to BIC: TIM2+F+G4\n",
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="IQ-TREE version 2.4.0")

    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: "iqtree2",
    )
    monkeypatch.setattr(
        phylogeny,
        "detect_iqtree",
        lambda: {"available": True, "version": "2.4.0", "reason": ""},
    )
    tree = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_IQTREE,
        iqtree_model_mode=phylogeny.IQTREE_MODEL_MFP,
        ufboot_replicates=1000,
        alrt_replicates=1000,
        run_fn=fake_run,
    )
    assert tree["status"] == "COMPUTED"
    assert tree["tool"] == "IQ-TREE"
    assert tree["selected_model"] == "TIM2+F+G4"
    assert tree["model_criterion"] == "BIC"
    methods = [item["method"] for item in tree["support"]["methods"]]
    assert phylogeny.SUPPORT_UFBOOT in methods
    assert phylogeny.SUPPORT_SH_ALRT in methods
    assert phylogeny.SUPPORT_FELSENSTEIN not in methods
    internals = [node for node in tree["nodes"] if not node["is_leaf"]]
    dual = [node for node in internals if node.get("support_by_method")]
    assert dual
    combo = dual[0]["support_by_method"]
    assert combo[phylogeny.SUPPORT_SH_ALRT] == pytest.approx(95.4)
    assert combo[phylogeny.SUPPORT_UFBOOT] == pytest.approx(98)
    assert dual[0]["support"] is None


def test_iqtree_malformed_treefile_is_parsing_error(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )

    def fake_run(argv, **kwargs):
        prefix = argv[argv.index("--prefix") + 1]
        Path(prefix + ".treefile").write_text("this is not newick", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="IQ-TREE version 2.4.0")

    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: "iqtree2",
    )
    monkeypatch.setattr(
        phylogeny,
        "detect_iqtree",
        lambda: {"available": True, "version": "2.4.0", "reason": ""},
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(
            result, method=phylogeny.METHOD_IQTREE, run_fn=fake_run
        )
    assert exc.value.category == "PARSING_ERROR"


def test_iqtree_timeout_category(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )

    def fake_run(argv, **kwargs):
        raise phylogeny.subprocess.TimeoutExpired(cmd=argv, timeout=1)

    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: "iqtree2",
    )
    monkeypatch.setattr(
        phylogeny,
        "detect_iqtree",
        lambda: {"available": True, "version": "2.4.0", "reason": ""},
    )
    with pytest.raises(phylogeny.PhylogenyError) as exc:
        phylogeny.infer_phylogeny(
            result, method=phylogeny.METHOD_IQTREE, run_fn=fake_run
        )
    assert exc.value.category == "TIMEOUT"


def test_iqtree_cache_stale_when_model_changes(monkeypatch):
    result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    validated = phylogeny.validate_msa_for_phylogeny(result)
    names = [item["tree_id"] for item in validated["leaf_records"]]
    newick = f"({names[0]}:0.1,{names[1]}:0.2,{names[2]}:0.3);"
    calls = {"n": 0}

    def fake_run(argv, **kwargs):
        calls["n"] += 1
        prefix = argv[argv.index("--prefix") + 1]
        Path(prefix + ".treefile").write_text(newick, encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="IQ-TREE version 2.4.0")

    monkeypatch.setattr(
        phylogeny.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: "iqtree2",
    )
    monkeypatch.setattr(
        phylogeny,
        "detect_iqtree",
        lambda: {"available": True, "version": "2.4.0", "reason": ""},
    )
    cache: dict = {}
    first = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_IQTREE,
        distance_model=phylogeny.DISTANCE_JC69,
        run_fn=fake_run,
        cache=cache,
    )
    second = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_IQTREE,
        distance_model=phylogeny.DISTANCE_JC69,
        run_fn=fake_run,
        cache=cache,
    )
    assert second["cache_status"] == "cached"
    third = phylogeny.infer_phylogeny(
        result,
        method=phylogeny.METHOD_IQTREE,
        iqtree_model_mode=phylogeny.IQTREE_MODEL_MFP,
        run_fn=fake_run,
        cache=cache,
    )
    assert third["cache_key"] != first["cache_key"]
    assert third["cache_status"] == "live"
    assert calls["n"] == 2


def test_dssp_maps_chain_and_skips_missing_residue():
    parsed = protein_structure.parse_dssp_legacy_output(
        (FIXTURES / "dssp_legacy_min.txt").read_text(encoding="utf-8")
    )
    assert parsed["residues"][0]["chain_id"] == "A"
    residues = [
        {"chain_id": "A", "auth_seq_id": 1, "one_letter": "I", "has_coordinates": True},
        {"chain_id": "A", "auth_seq_id": 2, "one_letter": "T", "has_coordinates": True},
        {"chain_id": "A", "auth_seq_id": 4, "one_letter": "G", "has_coordinates": True},
    ]
    mapped = protein_structure.map_dssp_assignments(parsed["residues"], residues)
    by_seq = {item["auth_seq_id"]: item for item in mapped["rows"]}
    assert by_seq[1]["mapped"] is True
    assert by_seq[1]["dssp_code"] == "H"
    assert by_seq[4]["mapped"] is False
    assert by_seq[4]["dssp_code"] is None
    assert mapped["n_unmapped_dssp"] == 1


def test_stride_injected_output_is_not_dssp():
    text = (
        "ASG  ILE A    1    1    H        360.00    360.00       45.00      0.00\n"
        "ASG  THR A    2    2    E        360.00    360.00       45.00      0.00\n"
    )
    result = protein_structure.assign_secondary_structure_stride(output_text=text)
    assert result["status"] == "EXPERIMENTAL"
    assert result["tool"] == "STRIDE"
    assert result["n_helix"] == 1
    assert result["n_sheet"] == 1
    missing = protein_structure.assign_secondary_structure_stride("")
    assert missing["status"] == "UNAVAILABLE"
    malformed = protein_structure.assign_secondary_structure_stride(
        output_text="returncode 0 but no ASG"
    )
    assert malformed["status"] == "ERROR"


def test_shrake_rupley_isolated_carbon_matches_4pi_r2():
    parsed = {
        "atoms": [
            {
                "group": "ATOM",
                "atom_name": "C",
                "element": "C",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "auth_asym_id": "A",
                "auth_seq_id": 1,
            }
        ]
    }
    result = protein_structure.shrake_rupley_sasa(parsed)
    radius = 1.70 + 1.4
    expected = 4.0 * math.pi * radius * radius
    assert result["status"] == "COMPUTED"
    assert result["sasa_angstrom2"] == pytest.approx(round(expected, 3), abs=0.001)
    assert result["not_msms"] is True


def test_surface_still_refuses_illustrative():
    with pytest.raises(protein_structure.StructureError) as exc:
        protein_structure.generate_molecular_surface(kind="illustrative")
    assert exc.value.category == "INVALID_INPUT"


def test_interchain_contacts_require_two_chains():
    parsed = {
        "atoms": [
            {
                "group": "ATOM",
                "atom_name": "CA",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
                "auth_asym_id": "A",
                "auth_seq_id": 1,
                "comp_id": "ALA",
            },
            {
                "group": "ATOM",
                "atom_name": "CA",
                "x": 4.0,
                "y": 0.0,
                "z": 0.0,
                "auth_asym_id": "B",
                "auth_seq_id": 1,
                "comp_id": "ALA",
            },
        ]
    }
    contacts = protein_structure.interchain_contacts(parsed, threshold_angstrom=8.0)
    assert contacts["status"] == "COMPUTED"
    assert contacts["n_pairs"] == 1
    same = protein_structure.interchain_contacts(
        {"atoms": [parsed["atoms"][0]]}, threshold_angstrom=8.0
    )
    assert same["n_pairs"] == 0
    assert same["status"] == "COMPUTED"


def test_ruleset2_and_deephf_never_return_zero():
    spacer = "GAGTCCGAGCAGAAGAAGAA"
    rs2 = crispr_ontarget.score_doench_ruleset2(spacer, nuclease="SpCas9")
    assert rs2["status"] == "UNAVAILABLE"
    assert rs2["score"] is None
    assert rs2["score"] != 0
    wrong = crispr_ontarget.score_doench_ruleset2(spacer, nuclease="AsCas12a")
    assert wrong["score"] is None
    short = crispr_ontarget.score_doench_ruleset2("ACGT", nuclease="SpCas9")
    assert short["score"] is None
    deephf = crispr_ontarget.score_deephf(spacer + "GGG", nuclease="SpCas9")
    assert deephf["score"] is None
    ranked = crispr.rank_guides(
        [{"guide_sequence": spacer, "pam_sequence": "GGG", "cas_system": "SpCas9"}]
    )
    assert ranked[0]["ruleset2_score"] is None
    assert ranked[0]["deephf_score"] is None
    assert ranked[0]["doench_score"] != ranked[0]["ruleset2_score"]


def test_tool_registry_does_not_scan_and_lists_iqtree():
    snap = tool_registry.collect_tool_snapshot(refresh=True)
    assert "IQ-TREE" in snap["tools"]
    assert "no recursive" in snap["scan"].lower()
    source = inspect.getsource(tool_registry.collect_tool_snapshot)
    assert "os.walk" not in source
    assert "C:\\\\" not in source


def test_local_blast_unavailable_without_executable(monkeypatch):
    monkeypatch.delenv("HELIXSCOPE_BLAST_DB", raising=False)
    monkeypatch.setattr(
        blast_search.tool_detection,
        "resolve_allowlisted_executable",
        lambda *args, **kwargs: None,
    )
    local = blast_search.local_blast_availability()
    assert local["available"] is False
    assert local["executables_detected"] is False
    database = local["database"]
    if database.get("available"):
        assert database.get("not_nt") is True
        assert database.get("not_nr") is True
        assert "nt/nr" in str(database.get("reason") or "").lower() or "declared" in str(
            database.get("reason") or ""
        ).lower()
    remote = blast_search.ncbi_blast_availability()
    assert remote["available"] is True
    assert remote["endpoint"] == blast_search.BLAST_ENDPOINT
    monkeypatch.setattr(blast_search, "declared_blast_prefix", lambda program="blastn": "")
    monkeypatch.setattr(
        blast_search,
        "_local_blast_database_record",
        lambda: {
            "available": False,
            "prefix": "",
            "reason": "HELIXSCOPE_BLAST_DB is unset and no tools/blast_db fixture is present.",
        },
    )
    hidden = blast_search.local_blast_availability()
    assert hidden["database"]["available"] is False
    assert hidden["available"] is False


def test_iqtree_runner_stays_shell_false():
    source = inspect.getsource(phylogeny._run_iqtree)
    assert "shell=False" in source
    assert "shell=True" not in source
    assert "-B" in source
    assert "-alrt" in source
