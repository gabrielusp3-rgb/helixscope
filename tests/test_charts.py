"""Testes pontuais de graficos Plotly contra os dados calculados."""

import pytest

from modules import crispr
from ui import charts


def test_guide_comparison_omits_missing_specificity():
    guides = crispr.evaluate_guides(
        crispr.find_guides("GAGTCCGAGCAGAAGAAGAAGGG", "SpCas9"),
        "GAGTCCGAGCAGAAGAAGAAGGG",
        "SpCas9",
        run_off_target=False,
    )
    figure = charts.guide_comparison_chart(guides)
    names = [trace.name for trace in figure.data]
    assert "Heuristic efficiency" in names
    assert "GC / 100" in names
    assert "Local specificity proxy" not in names
    assert "Specificity proxy" not in names
    efficiency = list(figure.data[0].y)
    assert efficiency[0] == guides[0]["doench_score"]


def test_entropy_and_kmer_charts_use_provided_values():
    profile = [{"midpoint": 10, "entropy": 1.5}, {"midpoint": 20, "entropy": 2.0}]
    entropy_fig = charts.entropy_profile_plot(profile)
    assert list(entropy_fig.data[0].x) == [10, 20]
    assert list(entropy_fig.data[0].y) == [1.5, 2.0]
    kmer_fig = charts.kmer_bar_chart(
        {"AAA": {"count": 3, "frequency": 50.0}, "TTT": {"count": 1, "frequency": 16.7}}
    )
    assert "AAA" in list(kmer_fig.data[0].y)
    with pytest.raises(ValueError):
        charts.entropy_profile_plot([])
    with pytest.raises(ValueError):
        charts.kmer_bar_chart({})
    with pytest.raises(ValueError):
        charts.sequence_feature_map(10, [{"name": "ORFs", "features": []}])


def test_gc_gauge_rejects_nan_and_out_of_range():
    figure = charts.gc_gauge(50.0)
    assert figure.data[0].value == pytest.approx(50.0)
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.gc_gauge(float("nan"))
    with pytest.raises(ValueError):
        charts.gc_gauge(150.0)
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.nucleotide_bar_chart({"A": {"count": 0, "frequency": 0.0}})


def test_alignment_column_map_uses_match_mismatch_gap():
    figure = charts.alignment_column_map(["match", "mismatch", "gap", "match"])
    names = {trace.name for trace in figure.data}
    assert names == {"match", "mismatch", "gap"}
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.alignment_column_map([])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.alignment_column_map(["match", "unknown"])


def test_blast_charts_require_real_hits():
    hits = [
        {
            "accession": "NM_000000",
            "bit_score": 40.1,
            "evalue": 0.00012,
            "identity_pct": 100.0,
            "query_coverage_pct": 100.0,
            "query_from": 1,
            "query_to": 20,
        }
    ]
    ranking = charts.blast_hit_ranking_chart(hits)
    assert list(ranking.data[0].x) == [40.1]
    evalue = charts.blast_evalue_chart(hits)
    assert list(evalue.data[0].y) == [0.00012]
    identity = charts.blast_identity_chart(hits)
    assert list(identity.data[0].y) == [100.0]
    coverage = charts.blast_coverage_chart(hits)
    assert list(coverage.data[0].y) == [100.0]
    overview = charts.blast_alignment_overview(hits, 20)
    assert list(overview.data[0].x) == [1, 20]
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.blast_hit_ranking_chart([])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.blast_evalue_chart([])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.blast_alignment_overview([], 20)


def test_msa_conservation_chart_uses_provided_scores():
    figure = charts.msa_conservation_chart([0, 1], [1.0, 0.5], ["A", "R"])
    assert list(figure.data[0].x) == [0, 1]
    assert list(figure.data[0].y) == [1.0, 0.5]
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.msa_conservation_chart([], [], [])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.msa_conservation_chart([0], [float("nan")], ["-"])


def test_phylogenetic_tree_figure_draws_real_edges_only():
    from tests.test_phylogeny import _msa_from_aligned
    from modules import phylogeny

    msa_result = _msa_from_aligned(
        ["a", "b", "c"],
        ["AAAAAAAAAA", "AAAAAAAATT", "TTTTTTTTTT"],
    )
    tree = phylogeny.infer_phylogeny(
        msa_result, method=phylogeny.METHOD_NJ, distance_model=phylogeny.DISTANCE_P
    )
    layout = tree["layouts"]["rectangular"]
    figure = charts.phylogenetic_tree_figure(
        layout,
        selected_tree_id=tree["leaves"][0]["tree_id"],
        rooted_label="unrooted",
    )
    line_traces = [
        trace for trace in figure.data if trace.mode == "lines" or trace.mode == "lines+text"
    ]
    topology_lines = [trace for trace in line_traces if trace.name != "Leaves"]
    assert len(topology_lines) >= len(layout["edges"])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.phylogenetic_tree_figure({"nodes": [], "edges": []})


def test_rna_structure_charts_use_provided_pairs():
    from modules import rna_folding

    sequence = "GGGAAACCC"
    pairs = [
        {
            "position_0based": 0,
            "base": "G",
            "paired_position_0based": 8,
            "paired_base": "C",
        }
    ]
    paths = rna_folding.arc_paths(9, pairs)
    figure = charts.rna_arc_diagram(sequence, paths, selected=0)
    assert figure.data[-1].text[0] == "G"
    assert figure.layout.dragmode == "zoom"
    layout = rna_folding.circular_layout(sequence, pairs)
    circular = charts.rna_circular_pairs(layout, selected=8)
    assert len(circular.data[-1].x) == 9
    assert circular.layout.dragmode == "zoom"
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.rna_arc_diagram("", [])
    with pytest.raises(ValueError, match="Insufficient data"):
        charts.rna_circular_pairs({"x": [], "y": [], "bases": [], "chords": []})


def test_dinucleotide_heatmap_does_not_plot_nan_as_zero():
    from modules import dna_analysis

    frequencies = dna_analysis.dinucleotide_frequencies("ATATATATAT")
    assert frequencies["CG"]["observed_expected"] != 0
    figure = charts.dinucleotide_heatmap(frequencies)
    z = list(figure.data[0].z)
    text = list(figure.data[0].text)
    cg = z[1][2]
    assert cg is None or (isinstance(cg, float) and cg != cg)
    assert "N/A" in {cell for row in text for cell in row}


def test_gc_skew_fill_does_not_treat_nan_as_zero_skew():
    figure = charts.gc_skew_plot([0.2, float("nan"), -0.1], window=10)
    fill_pos = list(figure.data[0].y)
    assert fill_pos[1] is None or (isinstance(fill_pos[1], float) and fill_pos[1] != fill_pos[1])


