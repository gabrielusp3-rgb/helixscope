"""Modos de RNA: folding global curto e analise local sem MFE global."""

from __future__ import annotations

import pytest

from modules import rna_folding


def test_folding_admission_keeps_global_and_local_distinct() -> None:
    short = rna_folding.folding_admission(80)
    assert short["global_fold"] == "AVAILABLE"
    assert short["local_long_rna"] == "AVAILABLE"
    long = rna_folding.folding_admission(150_000)
    assert long["global_fold"] == "RESOURCE_LIMIT"
    assert long["local_long_rna"] == "AVAILABLE"
    blocked = rna_folding.folding_admission(150_001)
    assert blocked["status"] == "RESOURCE_LIMIT"
    assert blocked["local_long_rna"] == "RESOURCE_LIMIT"


def test_global_fold_does_not_silently_shorten_a_long_rna() -> None:
    with pytest.raises(rna_folding.FoldingError) as caught:
        rna_folding.fold_rna("A" * (rna_folding.MAX_FOLD_NT + 1))
    assert caught.value.category == "RESOURCE_LIMIT"


@pytest.mark.skipif(
    not rna_folding.detect_python_bindings().get("available"),
    reason="ViennaRNA Python module RNA is not installed",
)
def test_local_analysis_uses_pfl_fold_on_the_full_short_rna() -> None:
    sequence = "GGGAAACCCUUUGGGAAACCC"
    result = rna_folding.local_long_rna_analysis(sequence, window=12, max_pair_span=10)
    assert result["mode"] == "LOCAL_LONG_RNA"
    assert result["global_structure"] is False
    assert result["input_length"] == len(sequence)
    assert result["analyzed_length"] == len(sequence)
    assert result["window_size"] == 12
    assert result["max_pair_span"] == 10
    assert result["engine"] == "ViennaRNA Python"
    assert result["engine_version"]
    assert result["pair_count"] >= 0
    assert result["visualized_length"] <= rna_folding.LOCAL_PREVIEW_PAIRS
    assert "global" in result["disclaimer"]
    if result["unpaired_probability_status"] == "COMPUTED":
        assert result["unpaired_probability_mean"] is not None
    else:
        assert result["unpaired_probability_mean"] is None


def test_ensemble_metrics_refuse_sequences_past_the_global_cap() -> None:
    with pytest.raises(rna_folding.FoldingError) as caught:
        rna_folding.ensemble_metrics("A" * (rna_folding.MAX_FOLD_NT + 1))
    assert caught.value.category == "RESOURCE_LIMIT"


@pytest.mark.skipif(
    not rna_folding.detect_python_bindings().get("available"),
    reason="ViennaRNA Python module RNA is not installed",
)
def test_ensemble_metrics_come_from_viennarna() -> None:
    result = rna_folding.ensemble_metrics("GGGAAACCC")
    assert result["status"] == "PREDICTED"
    assert result["mode"] == "GLOBAL_ENSEMBLE"
    assert result["analyzed_length"] == 9
    assert result["mfe_structure"]
    assert result["centroid_structure"]
    assert result["mea_structure"]
    assert result["ensemble_diversity"] is not None


def test_local_analysis_rejects_dna_and_overlong_input() -> None:
    with pytest.raises(rna_folding.FoldingError) as dna:
        rna_folding.local_long_rna_analysis("GGGAAACCCT")
    assert dna.value.category == "INVALID_INPUT"
    with pytest.raises(rna_folding.FoldingError) as over:
        rna_folding.local_long_rna_analysis("A" * (rna_folding.MAX_LOCAL_RNA_NT + 1))
    assert over.value.category == "RESOURCE_LIMIT"
