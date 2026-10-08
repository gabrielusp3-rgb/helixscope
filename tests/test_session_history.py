"""Historico de modulos e de analises. Nao grava disco e nao chama motores."""

from __future__ import annotations

from modules.session_history import (
    blank_history,
    clear_analyses,
    delete_analysis,
    go_back,
    go_forward,
    history_memory_chars,
    restore_analysis,
    sync_analysis_history,
    visit_module,
)

KNOWN = {"overview", "dna", "protein", "rna", "crispr", "msa", "phylogeny", "structure_3d", "compare"}


def test_two_sessions_do_not_share_or_delete_each_other() -> None:
    left = {"dna_input": "AAAAA"}
    right = {"dna_input": "CCCCC"}
    left_history = sync_analysis_history(left)
    right_history = sync_analysis_history(right)
    assert left_history is not right_history
    assert left_history["entries"][0]["input_hash"] != right_history["entries"][0]["input_hash"]
    delete_analysis(left_history, left_history["entries"][0]["id"])
    assert left_history["entries"] == []
    assert len(right_history["entries"]) == 1
    restored = restore_analysis(right_history, right_history["entries"][0]["id"], right)
    assert restored is not None
    assert right["dna_input"] == "CCCCC"

    back, forward, current = [], [], "overview"
    for target in ("dna", "protein", "rna"):
        back, forward, current = visit_module(back, forward, current, target, known=KNOWN)
    assert current == "rna"
    assert forward == []
    back, forward, current = go_back(back, forward, current, known=KNOWN)
    assert current == "protein"
    back, forward, current = go_back(back, forward, current, known=KNOWN)
    assert current == "dna"
    back, forward, current = go_back(back, forward, current, known=KNOWN)
    assert current == "overview"
    back, forward, current = go_back(back, forward, current, known=KNOWN)
    assert current == "overview"
    back, forward, current = go_forward(back, forward, current, known=KNOWN)
    assert current == "dna"


def test_forward_is_cleared_when_a_new_module_is_opened() -> None:
    back, forward, current = visit_module([], [], "overview", "dna", known=KNOWN)
    back, forward, current = visit_module(back, forward, current, "protein", known=KNOWN)
    back, forward, current = go_back(back, forward, current, known=KNOWN)
    assert current == "dna"
    assert forward == ["protein"]
    back, forward, current = visit_module(back, forward, current, "rna", known=KNOWN)
    assert current == "rna"
    assert forward == []


def test_same_module_does_not_grow_history() -> None:
    back, forward, current = visit_module(["overview"], [], "dna", "dna", known=KNOWN)
    assert back == ["overview"]
    assert current == "dna"


def test_analysis_history_delete_does_not_return_on_sync() -> None:
    state = {"dna_input": "ATGCATGC"}
    first = sync_analysis_history(state)
    assert len(first["entries"]) == 1
    entry_id = first["entries"][0]["id"]
    delete_analysis(first, entry_id)
    again = sync_analysis_history(state)
    assert again["entries"] == []


def test_clear_history_and_restore_retained_input() -> None:
    state = {"protein_input": "ACDEFGHIKL"}
    history = sync_analysis_history(state)
    entry_id = history["entries"][0]["id"]
    state["protein_input"] = "MMMM"
    restored = restore_analysis(history, entry_id, state)
    assert restored is not None
    assert state["protein_input"] == "ACDEFGHIKL"
    clear_analyses(history, ["protein|Protein profile|" + history["entries"][0]["input_hash"]])
    sync_analysis_history(state)
    assert history["entries"] == []


def test_history_keeps_one_hundred_short_entries_without_dropping_them() -> None:
    state: dict = {"helix_analysis_history": blank_history()}
    for index in range(100):
        state["dna_input"] = f"ACGT{index:04d}"
        state["helix_analysis_history"]["suppressed"] = []
        sync_analysis_history(state)
    history = state["helix_analysis_history"]
    assert len(history["entries"]) == 100
    used = history_memory_chars(history)
    state["dna_input"] = "TTTTTTTT"
    state["helix_analysis_history"]["suppressed"] = []
    sync_analysis_history(state)
    assert len(history["entries"]) == 100
    assert history["refused"]
    assert history_memory_chars(history) == used


def test_one_ten_and_fifty_entries_grow_without_a_jump() -> None:
    sizes = []
    for count in (1, 10, 50):
        state = {"helix_analysis_history": blank_history(), "rna_input": ""}
        for index in range(count):
            state["rna_input"] = "ACGU" * (index + 1)
            state["helix_analysis_history"]["suppressed"] = []
            sync_analysis_history(state)
        sizes.append(history_memory_chars(state["helix_analysis_history"]))
        assert len(state["helix_analysis_history"]["entries"]) == count
    assert sizes[0] < sizes[1] < sizes[2]
