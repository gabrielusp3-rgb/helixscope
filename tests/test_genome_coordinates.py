"""Coordinate roundtrip, PAM strand, display vs internal."""

from __future__ import annotations

import pytest

from modules import genome_coordinates


def test_internal_display_roundtrip_nonempty():
    for start, end in ((0, 1), (0, 23), (40, 63), (100, 250)):
        assert genome_coordinates.roundtrip_ok(start, end)
        shown = genome_coordinates.internal_to_display(start, end)
        back = genome_coordinates.display_to_internal(
            shown["start_1based"], shown["end_1based"]
        )
        assert back["start_0based"] == start
        assert back["end_0based"] == end


def test_empty_interval_has_no_display_bases():
    shown = genome_coordinates.internal_to_display(4, 4)
    assert shown["empty"] is True
    assert shown["start_1based"] is None


def test_display_rejects_zero_and_inverted():
    with pytest.raises(genome_coordinates.CoordinateError):
        genome_coordinates.display_to_internal(0, 10)
    with pytest.raises(genome_coordinates.CoordinateError):
        genome_coordinates.display_to_internal(10, 9)


def test_pam_plus_is_downstream_of_spacer():
    pam = genome_coordinates.pam_interval_on_strand(40, 60, strand="+", pam_nt=3)
    assert pam["pam_start_0based"] == 60
    assert pam["pam_end_0based"] == 63


def test_pam_minus_is_upstream_of_spacer_on_sense_contig():
    pam = genome_coordinates.pam_interval_on_strand(3, 23, strand="-", pam_nt=3)
    assert pam["pam_start_0based"] == 0
    assert pam["pam_end_0based"] == 3
    with pytest.raises(genome_coordinates.CoordinateError):
        genome_coordinates.pam_interval_on_strand(1, 21, strand="-", pam_nt=3)
