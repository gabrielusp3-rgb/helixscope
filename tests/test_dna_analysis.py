"""Testes unitarios do modulo modules.dna_analysis.

Cobrem validacao e deteccao de tipo, conteudo GC, temperatura de melting (Wallace
e salt-adjusted), complemento reverso com IUPAC, deteccao de ORFs e o tratamento
de erro do GC skew. Nenhum teste faz chamadas de rede.
"""

import math
from pathlib import Path

from tests.helix_apptest import open_module

import pytest

from modules import dna_analysis


@pytest.fixture
def dna_sequence() -> str:
    """Sequencia de DNA canonica reutilizavel nos testes."""
    return "ATGGCATTACGTACGTACGT"


@pytest.fixture
def gc_fifty_sequence() -> str:
    """Sequencia com conteudo GC conhecido de 50%."""
    return "GGCCAATT"


def test_validate_sequence_valid_dna():
    result = dna_analysis.validate_sequence("ACGTACGT")
    assert result["type"] == "DNA"
    assert result["is_valid"] is True
    assert result["invalid_chars"] == []


def test_validate_sequence_valid_rna():
    result = dna_analysis.validate_sequence("ACGUACGU")
    assert result["type"] == "RNA"
    assert result["is_valid"] is True


def test_validate_sequence_valid_protein():
    result = dna_analysis.validate_sequence("MKLPQERT")
    assert result["type"] == "PROTEIN"
    assert result["is_valid"] is True


def test_validate_sequence_invalid_chars():
    result = dna_analysis.validate_sequence("ACGT123")
    assert result["is_valid"] is False
    assert "1" in result["invalid_chars"]


def test_validate_for_molecule_dna_rejects_n():
    result = dna_analysis.validate_for_molecule("ACGTN", "DNA")
    assert result["is_valid"] is False
    assert "N" in result["invalid_chars"]


def test_validate_for_molecule_dna_rejects_rna_and_protein():
    rna = dna_analysis.validate_for_molecule("ACGUACGU", "DNA")
    assert rna["is_valid"] is False
    assert rna["type"] == "RNA"
    protein = dna_analysis.validate_for_molecule("MKLPQERT", "DNA")
    assert protein["is_valid"] is False
    assert protein["type"] == "PROTEIN"


def test_validate_for_molecule_rna_rejects_dna_and_protein():
    dna = dna_analysis.validate_for_molecule("ACGTACGT", "RNA")
    assert dna["is_valid"] is False
    assert dna["type"] == "DNA"
    assert "transcribe" in dna["rejection_reason"].lower() or "DNA" in dna["rejection_reason"]
    protein = dna_analysis.validate_for_molecule("MKLPQERT", "RNA")
    assert protein["is_valid"] is False
    assert protein["type"] == "PROTEIN"


def test_validate_for_molecule_protein_rejects_nucleic_acids():
    dna = dna_analysis.validate_for_molecule("ACGTACGT", "PROTEIN")
    assert dna["is_valid"] is False
    assert dna["type"] == "DNA"
    rna = dna_analysis.validate_for_molecule("ACGUACGU", "PROTEIN")
    assert rna["is_valid"] is False
    assert rna["type"] == "RNA"
    ok = dna_analysis.validate_for_molecule("MKLPQERT", "PROTEIN")
    assert ok["is_valid"] is True
    assert ok["type"] == "PROTEIN"


def test_validate_for_molecule_rejects_stop_and_x():
    result = dna_analysis.validate_for_molecule("MKVL*", "PROTEIN")
    assert result["is_valid"] is False
    assert "*" in result["invalid_chars"]
    unknown = dna_analysis.validate_for_molecule("MKVLX", "PROTEIN")
    assert unknown["is_valid"] is False
    assert "X" in unknown["invalid_chars"]


def test_gc_content_known_value(gc_fifty_sequence):
    expected = (4 / 8) * 100.0
    assert dna_analysis.gc_content(gc_fifty_sequence) == pytest.approx(expected)


def test_melting_temperature_wallace_short():
    expected = 2 * (1 + 1) + 4 * (1 + 1)
    assert dna_analysis.melting_temperature("ATGC") == pytest.approx(expected)


def test_melting_temperature_salt_adjusted_long():
    seq = "ACGTACGTACGTACGTACGT"
    na = 0.05
    gc_percent = 50.0
    n = len(seq)
    expected = 81.5 + 16.6 * math.log10(na) + 0.41 * gc_percent - 675.0 / n
    result = dna_analysis.melting_temperature(seq, na_concentration=na)
    assert result == pytest.approx(round(expected, 2), abs=0.01)


def test_reverse_complement_iupac():
    assert dna_analysis.reverse_complement("ATGC") == "GCAT"
    assert dna_analysis.reverse_complement("RYSWKM") == "KMWSRY"


def test_find_orfs_known_orf():
    orfs = dna_analysis.find_orfs("ATGAAATAA", min_length=9)
    assert any(
        orf["frame"] == 1
        and orf["start"] == 0
        and orf["end"] == 9
        and orf["protein"] == "MK"
        for orf in orfs
    )


def test_gc_skew_window_larger_than_sequence():
    with pytest.raises(ValueError):
        dna_analysis.gc_skew("ACGT", window=10)


@pytest.fixture
def nested_orfs() -> list:
    """ORFs com um par aninhado na fita direta e uma ORF na fita reversa."""
    return [
        {"frame": 1, "start": 0, "end": 300, "length_bp": 300, "protein": ""},
        {"frame": 1, "start": 90, "end": 300, "length_bp": 210, "protein": ""},
        {"frame": 2, "start": 400, "end": 600, "length_bp": 200, "protein": ""},
        {"frame": -1, "start": 0, "end": 300, "length_bp": 300, "protein": ""},
    ]


def test_select_non_overlapping_orfs_drops_nested(nested_orfs):
    kept = dna_analysis.select_non_overlapping_orfs(nested_orfs, min_length=150)
    spans = {(orf["frame"], orf["start"], orf["end"]) for orf in kept}
    assert (1, 0, 300) in spans
    assert (1, 90, 300) not in spans
    assert (2, 400, 600) in spans
    assert (-1, 0, 300) in spans


def test_select_non_overlapping_orfs_applies_min_length(nested_orfs):
    kept = dna_analysis.select_non_overlapping_orfs(nested_orfs, min_length=250)
    assert all(orf["length_bp"] >= 250 for orf in kept)
    assert len(kept) == 2


def test_select_non_overlapping_orfs_rejects_invalid_arguments(nested_orfs):
    with pytest.raises(ValueError):
        dna_analysis.select_non_overlapping_orfs(nested_orfs, min_length=0)
    with pytest.raises(ValueError):
        dna_analysis.select_non_overlapping_orfs(
            nested_orfs, max_overlap_fraction=1.5
        )


def test_nucleotide_composition_counts_all_symbols():
    composition = dna_analysis.nucleotide_composition("AACGTN")
    assert composition["A"]["count"] == 2
    assert composition["N"]["count"] == 1


def test_nucleotide_composition_counts_rna_u():
    composition = dna_analysis.nucleotide_composition("AAUU")
    assert composition["A"]["count"] == 2
    assert composition["U"]["count"] == 2
    assert composition["T"]["count"] == 0
    assert composition["A"]["frequency"] + composition["U"]["frequency"] == pytest.approx(100.0)


def test_molecular_weight_positive_for_dna():
    assert dna_analysis.molecular_weight("ATGC", "DNA") > 0.0


def test_find_restriction_sites_reports_known_enzyme():
    sites = dna_analysis.find_restriction_sites("AAAGAATTCAAA")
    assert sites["EcoRI"]["count"] == 1


def test_cpg_islands_returns_empty_for_short_sequence():
    assert dna_analysis.cpg_islands("ACGT", window=200) == []


def test_at_content_complements_gc_content(gc_fifty_sequence):
    at = dna_analysis.at_content(gc_fifty_sequence)
    gc = dna_analysis.gc_content(gc_fifty_sequence)
    assert at == pytest.approx(50.0)
    assert at + gc == pytest.approx(100.0)


def test_at_content_handles_rna_and_empty():
    assert dna_analysis.at_content("AAUU") == pytest.approx(100.0)
    assert dna_analysis.gc_content("AAUU") == pytest.approx(0.0)
    assert math.isnan(dna_analysis.at_content("NNNN"))
    assert math.isnan(dna_analysis.gc_content("NNNN"))
    assert math.isnan(dna_analysis.at_content(""))
    assert math.isnan(dna_analysis.gc_content(""))


def test_dinucleotide_frequencies_covers_all_pairs():
    frequencies = dna_analysis.dinucleotide_frequencies("ACGT")
    assert len(frequencies) == 16
    assert set(frequencies) == set(dna_analysis.DINUCLEOTIDES)
    assert sum(entry["count"] for entry in frequencies.values()) == 3


def test_dinucleotide_frequencies_cpg_odds_ratio():
    enriched = dna_analysis.dinucleotide_frequencies("CGCGCGCGCG")
    assert enriched["CG"]["count"] == 5
    assert enriched["CG"]["observed_expected"] == pytest.approx(2.222, abs=0.001)


def test_dinucleotide_frequencies_absent_base_returns_nan_odds():
    depleted = dna_analysis.dinucleotide_frequencies("ATATATATAT")
    assert depleted["CG"]["count"] == 0
    assert math.isnan(depleted["CG"]["observed_expected"])


def test_gc_sliding_window_known_profile():
    profile = dna_analysis.gc_sliding_window("GGGGCCCCAAAATTTT", window=4, step=4)
    assert [point["gc_percent"] for point in profile] == [100.0, 100.0, 0.0, 0.0]
    assert profile[0]["start"] == 0
    assert profile[0]["end"] == 4
    assert profile[1]["midpoint"] == 6


def test_gc_sliding_window_rejects_invalid_arguments():
    with pytest.raises(ValueError):
        dna_analysis.gc_sliding_window("ACGTACGT", window=0)
    with pytest.raises(ValueError):
        dna_analysis.gc_sliding_window("ACGT", window=100)


def test_melting_temperature_above_200_bp_returns_none():
    assert dna_analysis.melting_temperature("A" * 201) is None


def test_shannon_entropy_known_values():
    assert dna_analysis.shannon_entropy("AAAA") == pytest.approx(0.0)
    assert dna_analysis.shannon_entropy("ACGT") == pytest.approx(2.0)
    assert math.isnan(dna_analysis.shannon_entropy("NNNN"))
    assert math.isnan(dna_analysis.shannon_entropy(""))


def test_at_skew_known_values():
    assert dna_analysis.at_skew("AAAA") == pytest.approx(1.0)
    assert dna_analysis.at_skew("TTTT") == pytest.approx(-1.0)
    assert math.isnan(dna_analysis.at_skew("GGGG"))


def test_kmer_counts_overlapping_windows():
    counts = dna_analysis.kmer_counts("AAA", k=2)
    assert counts["AA"]["count"] == 2
    assert dna_analysis.kmer_counts("ACGT", k=1)["A"]["count"] == 1
    with pytest.raises(ValueError):
        dna_analysis.kmer_counts("ACGT", k=6)


def test_entropy_sliding_window_and_consistency():
    profile = dna_analysis.entropy_sliding_window("ACGTACGT", window=4, step=4)
    assert len(profile) == 2
    assert profile[0]["entropy"] == pytest.approx(2.0)
    assert profile[0]["midpoint"] == 3
    assert dna_analysis.composition_is_consistent("ACGTN")
    assert dna_analysis.entropy_sliding_window("ACGT", window=8, step=1) == []


def test_at_skew_windows_and_orf_coordinate_mapping():
    windows = dna_analysis.at_skew_windows("AAAATTTT", window=4)
    assert windows[0] == pytest.approx(1.0)
    assert windows[1] == pytest.approx(-1.0)
    mapped = dna_analysis.orf_coordinates_on_input(
        {"frame": 1, "start": 0, "end": 9}, 20
    )
    assert mapped == (0, 9)
    reverse_mapped = dna_analysis.orf_coordinates_on_input(
        {"frame": -1, "start": 0, "end": 9}, 20
    )
    assert reverse_mapped == (11, 20)
    with pytest.raises(ValueError):
        dna_analysis.orf_coordinates_on_input(
            {"frame": 1, "start": 0, "end": 50}, 20
        )


def test_parse_sequence_payload_single_fasta_and_raw():
    parsed = dna_analysis.parse_sequence_payload(">seq1\nATGCATGC\n")
    assert parsed["format"] == "fasta"
    assert parsed["record_count"] == 1
    assert parsed["sequence"] == "ATGCATGC"
    assert parsed["identifier"] == "seq1"
    raw = dna_analysis.parse_sequence_payload("  atgcatgc  ")
    assert raw["format"] == "raw"
    assert raw["sequence"] == "atgcatgc"
    empty = dna_analysis.parse_sequence_payload("   ")
    assert empty["format"] == "empty"
    assert empty["sequence"] == ""


def test_parse_sequence_payload_rejects_multifasta():
    payload = ">one\nAAAA\n>two\nCCCC\n"
    with pytest.raises(ValueError, match="2 records"):
        dna_analysis.parse_sequence_payload(payload)
    multi = dna_analysis.parse_fasta_records(payload)
    assert multi["record_count"] == 2
    assert [item["identifier"] for item in multi["records"]] == ["one", "two"]


def test_gc_sliding_window_all_n_is_unavailable_not_zero():
    profile = dna_analysis.gc_sliding_window("NNNN", window=2, step=2)
    assert len(profile) == 2
    assert all(math.isnan(point["gc_percent"]) for point in profile)


def test_melting_temperature_report_labels_methods():
    wallace = dna_analysis.melting_temperature_report("ATGC")
    assert wallace["method"] == "Wallace"
    assert wallace["status"] == "COMPUTED"
    assert wallace["value_c"] == pytest.approx(12.0)
    salt = dna_analysis.melting_temperature_report("ACGTACGTACGTACGTACGT")
    assert salt["method"] == "salt-adjusted"
    assert salt["status"] == "COMPUTED"
    long_seq = dna_analysis.melting_temperature_report("A" * 201)
    assert long_seq["status"] == "UNAVAILABLE"
    assert long_seq["value_c"] is None
    empty = dna_analysis.melting_temperature_report("NNNN")
    assert empty["status"] == "UNAVAILABLE"
    assert empty["value_c"] is None


def test_molecular_weight_nan_without_canonical_bases():
    assert math.isnan(dna_analysis.molecular_weight("NNNN", "DNA"))
    assert math.isnan(dna_analysis.molecular_weight("", "DNA"))


def test_gc_and_at_skew_denominator_zero_is_nan():
    assert math.isnan(dna_analysis.gc_skew_global("ATAT"))
    assert dna_analysis.gc_skew_global("GGGG") == pytest.approx(1.0)
    assert dna_analysis.gc_skew_global("CCCC") == pytest.approx(-1.0)
    assert math.isnan(dna_analysis.at_skew("GGCC"))


def test_windowed_profiles_coordinates_and_nan():
    profile = dna_analysis.windowed_profiles("GGGGCCCC", window=4, step=4)
    assert profile[0]["start"] == 0
    assert profile[0]["end"] == 4
    assert profile[0]["window_size"] == 4
    assert profile[0]["step_size"] == 4
    assert profile[0]["alphabet"] == "ACGT"
    assert profile[0]["gc_percent"] == pytest.approx(100.0)
    all_n = dna_analysis.windowed_profiles("NNNN", window=2, step=2)
    assert all(math.isnan(point["gc_percent"]) for point in all_n)
    assert all(math.isnan(point["gc_skew"]) for point in all_n)
    with pytest.raises(ValueError):
        dna_analysis.windowed_profiles("ACGT", window=0)


def test_kmer_summary_and_memory_cap():
    summary = dna_analysis.kmer_summary("AAA", k=2, top_n=5)
    assert summary["valid_windows"] == 2
    assert summary["unique_kmers"] == 1
    assert summary["top"][0]["kmer"] == "AA"
    assert summary["status"] == "COMPUTED"
    empty = dna_analysis.kmer_summary("NNN", k=2)
    assert empty["status"] == "UNAVAILABLE"
    assert math.isnan(empty["diversity"])
    with pytest.raises(ValueError, match="memoria|memory"):
        dna_analysis.kmer_counts("ACGT", k=1_000_000)
    with pytest.raises(ValueError):
        dna_analysis.kmer_counts("ACGT", k=dna_analysis.MAX_KMER_K + 1)


def test_supported_iupac_dna_symbols_and_reverse_complement():
    mapping = dna_analysis.supported_iupac_dna_symbols()
    assert set("ACGT") <= set(mapping)
    assert "U" not in mapping
    for symbol in mapping:
        rc = dna_analysis.reverse_complement(symbol)
        assert len(rc) == 1
        assert rc in dna_analysis.DNA_COMPLEMENT_IUPAC.values()
    assert dna_analysis.reverse_complement("U") == "A"
    with pytest.raises(ValueError):
        dna_analysis.reverse_complement("Z")


def test_find_orfs_predicted_status_and_input_coordinates():
    orfs = dna_analysis.find_orfs("ATGAAATAA", min_length=9)
    plus = [orf for orf in orfs if orf["frame"] == 1][0]
    assert plus["status"] == "PREDICTED"
    assert plus["strand"] == "+"
    assert plus["start_codon"] == "ATG"
    assert plus["stop_codon"] == "TAA"
    assert plus["input_start"] == 0
    assert plus["input_end"] == 9
    assert plus["start"] == 0
    assert plus["end"] == 9
    rc_source = dna_analysis.reverse_complement("ATGAAATAA")
    reverse_orfs = dna_analysis.find_orfs(rc_source, min_length=9)
    minus = [orf for orf in reverse_orfs if orf["strand"] == "-"]
    assert minus
    item = minus[0]
    assert item["input_start"] == 0
    assert item["input_end"] == 9
    mapped = dna_analysis.orf_coordinates_on_input(item, len(rc_source))
    assert mapped == (item["input_start"], item["input_end"])


def test_restriction_hits_include_enzyme_site_position_strand():
    sites = dna_analysis.find_restriction_sites("AAAGAATTCAAA")
    hits = sites["EcoRI"]["hits"]
    assert hits[0]["enzyme"] == "EcoRI"
    assert hits[0]["site"] == "GAATTC"
    assert hits[0]["position"] == 4
    assert hits[0]["strand"] == "+"


def test_cpg_dinucleotide_positions_are_in_sequence():
    seq = "AACGTTCG"
    positions = dna_analysis.cpg_dinucleotide_positions(seq)
    assert positions == [2, 6]
    for start in positions:
        assert 0 <= start < len(seq) - 1
        assert seq[start : start + 2] == "CG"
    assert dna_analysis.cpg_dinucleotide_positions("ATAT") == []


def test_aaaaa_hides_gc_window_when_stored_size_is_below_minimum():
    """Stored window 4 must not crash when the sequence becomes shorter than 10 bp."""
    app = open_module("dna", timeout=90)
    long_seq = "ACGT" * 20
    app.text_area(key="dna_text").set_value(long_seq).run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_gc_window").set_value(True).run()
    assert not app.exception
    app.session_state["dna_gc_window_size"] = 4
    app.text_area(key="dna_text").set_value("AAAAA").run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    labels = [str(getattr(item, "label", "")) for item in app.number_input]
    assert "Window size (bp)" not in labels
    body = " ".join(str(getattr(item, "value", "")) for item in app.markdown)
    assert "5 bp" in body or "5" in body
    app.text_area(key="dna_text").set_value(long_seq).run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_gc_window").set_value(True).run()
    assert not app.exception
    labels = [str(getattr(item, "label", "")) for item in app.number_input]
    assert "Window size (bp)" in labels


def test_short_dna_clamps_a_stale_sequence_viewer_index():
    """A jump index from a long sequence must not exceed the new length."""
    app = open_module("dna", timeout=90)
    long_seq = "ACGT" * 20
    app.text_area(key="dna_text").set_value(long_seq).run()
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    app.toggle(key="dna_seq_view").set_value(True).run()
    assert not app.exception
    app.session_state["dna_view_jump"] = 400
    app.text_area(key="dna_text").set_value("AAAAA").run()
    app.button(key="dna_analyze").click().run()
    app.toggle(key="dna_seq_view").set_value(True).run()
    assert not app.exception
    jump = app.number_input(key="dna_view_jump")
    assert int(jump.value) <= 4


def test_dna_analyze_does_not_crash_on_stale_result_guard():
    app = open_module("dna", timeout=60)
    assert not app.exception
    app.text_area(key="dna_text").set_value("ATGGCATTACGTACGTACGT").run()
    assert not app.exception
    app.button(key="dna_analyze").click().run()
    assert not app.exception
    errors = " ".join(str(getattr(item, "value", item)) for item in list(app.error) or [])
    assert "hashes_match" not in errors
    assert "AttributeError" not in errors
