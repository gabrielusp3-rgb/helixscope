"""Metadados taxonomicos reais anexados a folhas de uma arvore.

A taxonomia nao altera a topologia. Falha de recuperacao NCBI nao destroi a
inferencia filogenetica. Organismo declarado pelo usuario nunca e inventado.
Uma especie pode ter varias sequencias; folhas nao sao fundidas.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import hashlib
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence

from . import ncbi_fetch, provenance

TAXONOMY_SOURCE_NCBI: str = "NCBI Taxonomy"
TAXONOMY_SOURCE_DECLARED: str = "declared_on_record"
TAXONOMY_SOURCE_USER: str = "user-provided"


class TaxonomyError(Exception):
    """Falha classificada de taxonomia, isolada da filogenia.

    Attributes:
        category: INVALID_INPUT, NOT_FOUND, SOURCE_UNAVAILABLE, UNAVAILABLE.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "UNAVAILABLE")


def cache_key(query: str, source: str = TAXONOMY_SOURCE_NCBI) -> str:
    """Chave de cache da taxonomia, separada da chave da arvore.

    Args:
        query: Taxon ID ou nome cientifico.
        source: Fonte (NCBI Taxonomy).

    Returns:
        Hex SHA-256.

    Raises:
        Nenhum.
    """
    payload = "|".join(
        [str(source), str(query).strip().lower(), provenance.HELIXSCOPE_VERSION]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def declared_from_msa_row(row: Mapping[str, Any]) -> dict:
    """Metadados de organismo ja presentes na linha do MSA, sem fetch.

    Args:
        row: Linha do envelope MSA.

    Returns:
        Dict com organism_declared, accession, version, source. taxonomy e
        None. Nomes comuns nao sao inventados.

    Raises:
        Nenhum.
    """
    organism = str(row.get("organism") or "").strip()
    source = str(row.get("source") or "").strip()
    return {
        "organism_declared": organism,
        "organism_display": organism if organism else "unknown / user-provided",
        "accession": str(row.get("accession") or ""),
        "version": str(row.get("version") or ""),
        "sequence_hash": str(row.get("hash") or ""),
        "identifier": str(row.get("identifier") or ""),
        "metadata_source": source or TAXONOMY_SOURCE_USER,
        "taxonomy": None,
        "status": "DECLARED" if organism else "UNAVAILABLE",
    }


def attach_declared_metadata(tree_result: Mapping[str, Any]) -> dict:
    """Copia o envelope da arvore e preenche organismo declarado por folha.

    Args:
        tree_result: Envelope de phylogeny.infer_phylogeny.

    Returns:
        Novo dict. Topologia, Newick e tree_hash inalterados.

    Raises:
        TaxonomyError: INVALID_INPUT se nao houver folhas.
    """
    leaves = list(tree_result.get("leaves") or [])
    if not leaves:
        raise TaxonomyError("Tree has no leaves to annotate.", "INVALID_INPUT")
    updated = dict(tree_result)
    updated["leaves"] = [dict(leaf) for leaf in leaves]
    updated["taxonomy_layer"] = {
        "kind": "taxonomic_data",
        "not_phylogenetic_tree": True,
        "status": "DECLARED",
        "source": TAXONOMY_SOURCE_DECLARED,
        "note": (
            "Declared organism/accession from the MSA members. "
            "This is not a taxonomic tree and does not change phylogeny."
        ),
    }
    return updated


def attach_ncbi_taxonomy(
    tree_result: Mapping[str, Any],
    *,
    email: str,
    api_key: Optional[str] = None,
    cache: Optional[Dict[str, dict]] = None,
    fetch_fn: Optional[Callable[..., Mapping[str, Any]]] = None,
) -> dict:
    """Busca taxonomia NCBI para organismos declarados. Falhas ficam por folha.

    Args:
        tree_result: Envelope da arvore (topologia ja calculada).
        email: E-mail exigido pelo Entrez.
        api_key: Chave NCBI opcional.
        cache: Cache separado da arvore, chaveado por cache_key.
        fetch_fn: ncbi_fetch.fetch_taxonomy injetavel (testes).

    Returns:
        Novo envelope. tree_hash e Newick nao mudam. Folhas sem organismo
        permanecem unknown / user-provided. Erro NCBI nao levanta para fora
        se ao menos a arvore existir: cada folha registra o erro.

    Raises:
        TaxonomyError: INVALID_INPUT se email vazio ou arvore sem folhas.
            Nao converte falha de uma folha em destruição da arvore.

    Nota biologica:
        Taxon ID, nome cientifico, linhagem e rank vem do banco Taxonomy.
        Nao e a arvore filogenetica calculada. Mesma especie, varias folhas.
    """
    if not str(email or "").strip():
        raise TaxonomyError("NCBI Taxonomy requires a contact email.", "INVALID_INPUT")
    leaves = list(tree_result.get("leaves") or [])
    if not leaves:
        raise TaxonomyError("Tree has no leaves to annotate.", "INVALID_INPUT")
    getter = fetch_fn or ncbi_fetch.fetch_taxonomy
    store = cache if cache is not None else {}
    annotated: List[dict] = []
    failures: List[dict] = []
    retrieved = 0
    for leaf in leaves:
        item = dict(leaf)
        query = str(item.get("organism_declared") or "").strip()
        if not query:
            item["taxonomy"] = None
            item["organism_display"] = "unknown / user-provided"
            annotated.append(item)
            continue
        key = cache_key(query)
        record: Optional[Mapping[str, Any]] = store.get(key)
        if record is None:
            try:
                record = getter(query, email, api_key=api_key)
                store[key] = dict(record)
            except (ncbi_fetch.NCBIQueryError, ValueError, RuntimeError) as exc:
                failures.append(
                    {
                        "tree_id": item.get("tree_id"),
                        "query": query,
                        "error": str(exc),
                    }
                )
                item["taxonomy"] = {
                    "status": "UNAVAILABLE",
                    "source": TAXONOMY_SOURCE_NCBI,
                    "taxon_id": "",
                    "scientific_name": "",
                    "rank": "",
                    "lineage": [],
                    "common_name": "",
                    "reason": str(exc),
                }
                annotated.append(item)
                continue
        item["taxonomy"] = _public_taxonomy(record)
        scientific = str((item["taxonomy"] or {}).get("scientific_name") or "")
        if scientific:
            item["organism_display"] = scientific
        retrieved += 1
        annotated.append(item)
    updated = dict(tree_result)
    updated["leaves"] = annotated
    updated["taxonomy_layer"] = {
        "kind": "taxonomic_data",
        "not_phylogenetic_tree": True,
        "status": "RETRIEVED" if retrieved else "UNAVAILABLE",
        "source": TAXONOMY_SOURCE_NCBI,
        "n_retrieved": retrieved,
        "n_failed": len(failures),
        "failures": failures,
        "tree_hash_unchanged": updated.get("tree_hash"),
        "note": (
            "NCBI Taxonomy annotations. Classification, not phylogenetic inference. "
            "Updating taxonomy does not recompute the tree."
        ),
    }
    return updated


def color_groups_by_rank(
    tree_result: Mapping[str, Any],
    rank: str = "species",
) -> dict:
    """Agrupa folhas por rank taxonomico real para coloracao visual.

    Args:
        tree_result: Envelope com taxonomy por folha.
        rank: species, genus, family, order, class, phylum.

    Returns:
        Dict leaf_tree_id -> grupo, legend (grupo -> n). Ausencia = Unknown.

    Raises:
        TaxonomyError: INVALID_INPUT se o rank nao for suportado.

    Nota biologica:
        Cor e visual. Nao e distancia filogenetica nem conservacao de MSA.
    """
    allowed = {"species", "genus", "family", "order", "class", "phylum"}
    chosen = str(rank or "species").strip().lower()
    if chosen not in allowed:
        raise TaxonomyError(
            "Taxonomy coloring rank must be species, genus, family, order, class or phylum.",
            "INVALID_INPUT",
        )
    groups: Dict[str, str] = {}
    counts: Dict[str, int] = {}
    for leaf in list(tree_result.get("leaves") or []):
        label = _rank_label(leaf, chosen)
        groups[str(leaf.get("tree_id") or "")] = label
        counts[label] = counts.get(label, 0) + 1
    return {
        "rank": chosen,
        "groups": groups,
        "legend": counts,
        "disclaimer": (
            "Leaf colors follow NCBI taxonomic rank when retrieved. "
            "Color is not phylogenetic distance and not MSA conservation."
        ),
        "status": "COMPUTED",
    }


def _public_taxonomy(record: Mapping[str, Any]) -> dict:
    common = str(record.get("common_name") or "").strip()
    return {
        "status": "RETRIEVED",
        "source": str(record.get("source") or TAXONOMY_SOURCE_NCBI),
        "taxon_id": str(record.get("taxon_id") or ""),
        "scientific_name": str(record.get("scientific_name") or ""),
        "rank": str(record.get("rank") or ""),
        "lineage": list(record.get("lineage") or []),
        "common_name": common,
        "retrieved_at_utc": str(record.get("retrieved_at_utc") or ""),
        "reason": "",
    }


def _rank_label(leaf: Mapping[str, Any], rank: str) -> str:
    taxonomy = leaf.get("taxonomy") or {}
    if not taxonomy or str(taxonomy.get("status") or "") != "RETRIEVED":
        return "Unknown"
    if rank == "species":
        name = str(taxonomy.get("scientific_name") or "").strip()
        return name or "Unknown"
    for item in list(taxonomy.get("lineage") or []):
        if str(item.get("rank") or "").lower() == rank:
            name = str(item.get("scientific_name") or "").strip()
            if name:
                return name
    current_rank = str(taxonomy.get("rank") or "").lower()
    if current_rank == rank:
        return str(taxonomy.get("scientific_name") or "").strip() or "Unknown"
    return "Unknown"
