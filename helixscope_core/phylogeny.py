"""Phylogeny facade. Does not add methods. NJ/UPGMA/IQ-TREE/FastTree only."""

from __future__ import annotations

from typing import Any, Mapping

from modules import phylogeny
from modules.phylogeny import (
    DEFAULT_BOOTSTRAP_SEED,
    DISTANCE_JC69,
    DISTANCE_P,
    METHOD_FASTTREE,
    METHOD_IQTREE,
    METHOD_NJ,
    METHOD_UPGMA,
    infer_phylogeny,
    parse_newick,
)

__all__ = (
    "DEFAULT_BOOTSTRAP_SEED",
    "DISTANCE_JC69",
    "DISTANCE_P",
    "METHOD_FASTTREE",
    "METHOD_IQTREE",
    "METHOD_NJ",
    "METHOD_UPGMA",
    "infer_phylogeny",
    "infer_tree",
    "parse_newick",
    "annotate_taxonomy",
)


def infer_tree(msa_result: Mapping[str, Any], *, method: str = METHOD_NJ, **kwargs: Any) -> dict:
    """Infer a tree from a completed MSA. Keyword-only engine parameters."""
    return infer_phylogeny(msa_result, method=method, **kwargs)


def annotate_taxonomy(
    tree_result: Mapping[str, Any],
    *,
    email: str,
    api_key: str | None = None,
    declared_organisms: Mapping[str, str] | None = None,
) -> dict:
    """Attach NCBI taxonomy to leaves. Does not recompute topology.

    Args:
        tree_result: Completed phylogeny envelope.
        email: Entrez contact.
        api_key: Optional NCBI key.
        declared_organisms: Optional tree_id/leaf id -> organism name.

    Returns:
        New envelope. Newick and tree_hash must remain unchanged.

    Raises:
        TaxonomyError: INVALID_INPUT if email empty or tree has no leaves.
    """
    from modules import taxonomy

    updated = dict(tree_result)
    leaves = [dict(leaf) for leaf in list(updated.get("leaves") or [])]
    declared = dict(declared_organisms or {})
    if declared:
        for leaf in leaves:
            key = str(leaf.get("tree_id") or leaf.get("id") or leaf.get("identifier") or "")
            name = str(declared.get(key) or "").strip()
            if name:
                leaf["organism_declared"] = name
        updated["leaves"] = leaves
    annotated = taxonomy.attach_ncbi_taxonomy(updated, email=email, api_key=api_key)
    original_hash = tree_result.get("tree_hash")
    original_newick = tree_result.get("newick")
    if original_hash is not None:
        annotated["tree_hash"] = original_hash
        layer = dict(annotated.get("taxonomy_layer") or {})
        layer["tree_hash_unchanged"] = original_hash
        annotated["taxonomy_layer"] = layer
    if original_newick is not None:
        annotated["newick"] = original_newick
    annotated["taxonomy_is_not_phylogeny"] = True
    return annotated
