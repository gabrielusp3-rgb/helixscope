"""Complexos macromoleculares depositados (Cas-guia-DNA) sem geometria sintetica.

Classifica chains a partir de polymer_type e da descricao depositada no mmCIF.
Nao monta Cas9+guide+DNA procedural. Coordenadas so as do ficheiro RCSB.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence

from . import protein_structure, provenance

ROLE_CAS_PROTEIN: str = "cas_protein"
ROLE_GUIDE_RNA: str = "guide_rna"
ROLE_TARGET_DNA: str = "target_dna"
ROLE_NONTARGET_DNA: str = "nontarget_dna"
ROLE_DNA: str = "dna"
ROLE_RNA: str = "rna"
ROLE_PROTEIN: str = "protein"
ROLE_UNKNOWN: str = "unknown"

GUIDE_TOKENS: tuple[str, ...] = ("sgrna", "guide rna", "guide-rna", "crRNA", "tracr")
TARGET_TOKENS: tuple[str, ...] = ("target dna", "target strand", "target dna strand")
NONTARGET_TOKENS: tuple[str, ...] = (
    "non-target",
    "nontarget",
    "non target",
    "pam-containing",
)
CAS_TOKENS: tuple[str, ...] = ("cas9", "cas 9", "cas12", "csn1", "endonuclease cas")


class ComplexError(Exception):
    """Falha classificada de complexo estrutural.

    Attributes:
        category: INVALID_INPUT, PARSING_ERROR, UNAVAILABLE, RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str) -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT")


def classify_chain_role(chain: Mapping[str, Any], description: str = "") -> str:
    """Atribui um papel a partir do tipo de polimero e da descricao depositada.

    Args:
        chain: Chain parseada (polymer_type, sequence).
        description: _entity.pdbx_description quando existir.

    Returns:
        Papel. Sem descricao, DNA/RNA/proteina genericos — nao inventa
        target vs non-target.

    Raises:
        Nenhum.

    Nota biologica:
        Os tokens vêm de descricoes PDB reais (ex.: 4UN3 SGRNA, TARGET DNA
        STRAND, NON-TARGET DNA STRAND, Cas9). Nao sao uma inferencia
        estrutural propria.
    """
    ptype = str(chain.get("polymer_type") or "").lower()
    blob = f"{description} {chain.get('description') or ''}".lower()
    if any(token.lower() in blob for token in CAS_TOKENS):
        return ROLE_CAS_PROTEIN
    if ptype.startswith("polypeptide"):
        return ROLE_PROTEIN
    if "polyribonucleotide" in ptype or ptype == "rna":
        if any(token.lower() in blob for token in GUIDE_TOKENS):
            return ROLE_GUIDE_RNA
        return ROLE_RNA
    if "polydeoxyribonucleotide" in ptype or ptype == "dna":
        if any(token.lower() in blob for token in NONTARGET_TOKENS):
            return ROLE_NONTARGET_DNA
        if any(token.lower() in blob for token in TARGET_TOKENS):
            return ROLE_TARGET_DNA
        return ROLE_DNA
    return ROLE_UNKNOWN


def map_guide_to_rna_chain(
    guide_sequence: str,
    chains: Sequence[Mapping[str, Any]],
) -> dict:
    """Procura o protoespacador do guia na chain de RNA depositada.

    Args:
        guide_sequence: 20 nt do guia (sem PAM).
        chains: Chains ja classificadas.

    Returns:
        Dict status (EXACT, UNMAPPED, UNAVAILABLE), chain_id, start_0based.
        Sem match: UNMAPPED, nao um highlight inventado.

    Raises:
        ComplexError: INVALID_INPUT se o guia for vazio.

    Nota biologica:
        O complexo depositado tem o guia da experiencia. O guia do usuario
        so mapeia se a subsequencia existir na chain de RNA.
    """
    spacer = str(guide_sequence or "").strip().upper().replace("T", "U")
    if not spacer:
        raise ComplexError("Guide sequence is empty.", "INVALID_INPUT")
    rna_chains = [
        item
        for item in chains
        if str(item.get("role") or "") in {ROLE_GUIDE_RNA, ROLE_RNA}
    ]
    if not rna_chains:
        return {
            "status": "UNAVAILABLE",
            "reason": "This structure has no RNA chain to map a guide.",
            "chain_id": "",
            "start_0based": None,
            "end_0based": None,
        }
    for chain in rna_chains:
        polymer = str(chain.get("sequence") or "").upper().replace("T", "U")
        index = polymer.find(spacer)
        if index >= 0:
            return {
                "status": "EXACT",
                "chain_id": str(chain.get("chain_id") or ""),
                "start_0based": index,
                "end_0based": index + len(spacer),
                "method": "Exact subsequence match on deposited RNA polymer sequence",
            }
    return {
        "status": "UNMAPPED",
        "reason": (
            "The current guide spacer is not an exact subsequence of the "
            "deposited RNA chain. HelixScope does not force a highlight."
        ),
        "chain_id": str(rna_chains[0].get("chain_id") or ""),
        "start_0based": None,
        "end_0based": None,
    }


def map_target_to_dna_chain(
    target_sequence: str,
    chains: Sequence[Mapping[str, Any]],
) -> dict:
    """Procura o DNA alvo por identidade de sequencia, nao por ordem de chain.

    Args:
        target_sequence: DNA alvo do usuario.
        chains: Chains ja classificadas.

    Returns:
        EXACT, UNMAPPED ou UNAVAILABLE.

    Raises:
        ComplexError: INVALID_INPUT se o alvo for vazio.

    Nota biologica:
        Uma chain de DNA presente nao e automaticamente o alvo. Sem match
        de sequencia o status e UNMAPPED.
    """
    needle = str(target_sequence or "").strip().upper().replace("U", "T")
    if not needle:
        raise ComplexError("Target DNA sequence is empty.", "INVALID_INPUT")
    dna_chains = [
        item
        for item in chains
        if str(item.get("role") or "") in {ROLE_TARGET_DNA, ROLE_DNA, ROLE_NONTARGET_DNA}
    ]
    if not dna_chains:
        return {
            "status": "UNAVAILABLE",
            "reason": "This structure has no DNA chain to map a target.",
            "chain_id": "",
            "start_0based": None,
            "end_0based": None,
        }
    preferred = [
        item for item in dna_chains if str(item.get("role") or "") == ROLE_TARGET_DNA
    ] or dna_chains
    for chain in preferred:
        polymer = str(chain.get("sequence") or "").upper().replace("U", "T")
        index = polymer.find(needle)
        if index >= 0:
            return {
                "status": "EXACT",
                "chain_id": str(chain.get("chain_id") or ""),
                "start_0based": index,
                "end_0based": index + len(needle),
                "method": "Exact subsequence match on deposited DNA polymer sequence",
            }
    return {
        "status": "UNMAPPED",
        "reason": (
            "The current target DNA is not an exact subsequence of the "
            "deposited DNA chain. HelixScope does not force a highlight."
        ),
        "chain_id": str(preferred[0].get("chain_id") or ""),
        "start_0based": None,
        "end_0based": None,
    }


def map_pam_on_nontarget(
    pam: str,
    chains: Sequence[Mapping[str, Any]],
) -> dict:
    """Localiza o PAM na chain de DNA non-target se a sequencia depositada o contiver.

    Args:
        pam: Motivo (ex. NGG nao se expande; usa GG se o PAM for NGG).
        chains: Chains classificadas.

    Returns:
        Dict status, chain_id, indices. NGG: procura GG na chain nontarget/DNA.

    Raises:
        Nenhum.
    """
    motif = str(pam or "").strip().upper().replace("U", "T")
    if motif in {"NGG", ""}:
        needle = "GG"
        method = "Exact GG search on deposited non-target/DNA chain (SpCas9 PAM NGG)"
    else:
        if "N" in motif:
            return {
                "status": "UNAVAILABLE",
                "reason": "Ambiguous PAM letters other than NGG are not expanded into fake sites.",
                "chain_id": "",
                "start_0based": None,
                "end_0based": None,
            }
        needle = motif
        method = "Exact PAM string on deposited DNA chain"
    dna_chains = [
        item
        for item in chains
        if str(item.get("role") or "") in {ROLE_NONTARGET_DNA, ROLE_DNA, ROLE_TARGET_DNA}
    ]
    preferred = [
        item for item in dna_chains if str(item.get("role") or "") == ROLE_NONTARGET_DNA
    ] or dna_chains
    if not preferred:
        return {
            "status": "UNAVAILABLE",
            "reason": "This structure has no DNA chain to map a PAM.",
            "chain_id": "",
            "start_0based": None,
            "end_0based": None,
        }
    for chain in preferred:
        polymer = str(chain.get("sequence") or "").upper().replace("U", "T")
        index = polymer.find(needle)
        if index >= 0:
            return {
                "status": "EXACT",
                "chain_id": str(chain.get("chain_id") or ""),
                "start_0based": index,
                "end_0based": index + len(needle),
                "method": method,
            }
    return {
        "status": "UNMAPPED",
        "reason": "The PAM motif is not an exact subsequence of the deposited DNA chain.",
        "chain_id": str(preferred[0].get("chain_id") or ""),
        "start_0based": None,
        "end_0based": None,
    }


def annotate_parsed_complex(
    parsed: Mapping[str, Any],
    *,
    structure_id: str,
    source: str,
    metadata: Optional[Mapping[str, Any]] = None,
    guide_sequence: str = "",
    pam: str = "",
) -> dict:
    """Empacota um mmCIF/PDB ja parseado como complexo. Sem XYZ novos.

    Args:
        parsed: Saida de protein_structure.parse_structure_text.
        structure_id: PDB ID.
        source: RCSB PDB.
        metadata: Entry metadata.
        guide_sequence: Guia do usuario, opcional.
        pam: PAM do usuario, opcional.

    Returns:
        Envelope COMPLEX com chains, roles, mappings de guia/PAM.

    Raises:
        ComplexError: PARSING_ERROR se nao houver chains.
    """
    chains = [dict(item) for item in list(parsed.get("chains") or [])]
    if not chains:
        raise ComplexError("Structure has no polymer chains.", "PARSING_ERROR")
    descriptions = dict(parsed.get("entity_descriptions") or {})
    annotated: List[dict] = []
    for chain in chains:
        entity_id = str(chain.get("entity_id") or "")
        description = str(
            chain.get("description")
            or descriptions.get(entity_id)
            or descriptions.get(str(chain.get("chain_id") or ""))
            or ""
        )
        role = classify_chain_role(chain, description)
        item = dict(chain)
        item["role"] = role
        item["description"] = description
        item["molecule"] = _molecule_from_role(role, str(chain.get("polymer_type") or ""))
        annotated.append(item)
    guide_map = (
        map_guide_to_rna_chain(guide_sequence, annotated)
        if guide_sequence
        else {"status": "UNAVAILABLE", "reason": "No user guide supplied.", "chain_id": ""}
    )
    pam_map = (
        map_pam_on_nontarget(pam, annotated)
        if pam or guide_sequence
        else {"status": "UNAVAILABLE", "reason": "No PAM supplied.", "chain_id": ""}
    )
    meta = dict(metadata or {})
    kind = str(meta.get("kind") or parsed.get("kind") or "experimental")
    envelope = provenance.analysis_envelope(
        module="COMPLEX_STRUCTURE",
        payload={
            "structure_id": str(structure_id).upper(),
            "n_chains": len(annotated),
            "kind": kind,
        },
        status="RETRIEVED" if kind == "experimental" else "PREDICTED",
        algorithm="Deposited macromolecular complex (mmCIF/PDB parser)",
        parameters={"structure_id": str(structure_id).upper()},
        source=source,
        accession=str(structure_id).upper(),
    )
    envelope.update(
        {
            "status": "RETRIEVED" if kind == "experimental" else "PREDICTED",
            "kind": kind,
            "kind_label": protein_structure.structure_kind_label(kind),
            "structure_id": str(structure_id).upper(),
            "source": source,
            "method": str(meta.get("method") or parsed.get("method") or ""),
            "resolution_angstrom": meta.get("resolution_angstrom", parsed.get("resolution_angstrom")),
            "chains": annotated,
            "atoms": list(parsed.get("atoms") or []),
            "n_atoms": len(list(parsed.get("atoms") or [])),
            "content_hash": parsed.get("content_hash") or "",
            "guide_mapping": guide_map,
            "pam_mapping": pam_map,
            "target_mapping": {"status": "UNAVAILABLE", "reason": "Use map_target_to_dna_chain with an explicit target sequence."},
            "sequence_identities": {
                item.get("role"): provenance.sequence_digest(str(item.get("sequence") or ""))
                for item in annotated
                if item.get("sequence")
            },
            "mapping_status": {
                "guide": guide_map.get("status"),
                "pam": pam_map.get("status"),
            },
            "retrieved_at": provenance.utc_now(),
            "synthetic_complex": False,
            "disclaimer": (
                "Experimental/predicted coordinates from the deposited file. "
                "HelixScope did not assemble Cas9, guide and DNA procedurally. "
                "Guide and PAM highlights require an exact sequence match on "
                "the deposited polymer."
            ),
        }
    )
    return envelope


def annotate_from_polymer_entities(
    entities: Sequence[Mapping[str, Any]],
    *,
    structure_id: str,
    metadata: Optional[Mapping[str, Any]] = None,
    guide_sequence: str = "",
    target_sequence: str = "",
    pam: str = "",
) -> dict:
    """Classifica cadeias a partir do JSON polymer_entity da Data API.

    Args:
        entities: Saidas de parse_rcsb_polymer_entity (sem XYZ).
        structure_id: PDB ID.
        metadata: Entry metadata (metodo, resolucao).
        guide_sequence: Guia do usuario.
        target_sequence: DNA alvo do usuario.
        pam: PAM.

    Returns:
        Envelope RETRIEVED de identidade. Coordinates UNAVAILABLE ate o mmCIF.

    Raises:
        ComplexError: PARSING_ERROR se vazio.

    Nota biologica:
        Papéis vêm da descricao e do tipo de polimero depositados (4UN3:
        SGRNA, Cas9/Csn1, TARGET DNA STRAND, NON-TARGET DNA STRAND), nao da
        ordem alfabetica do chain ID.
    """
    chains: List[dict] = []
    for entity in entities:
        parsed = protein_structure.parse_rcsb_polymer_entity(entity) if "entity_poly" in entity else dict(entity)
        chain_id = ""
        auth_ids = list(parsed.get("auth_asym_ids") or [])
        if auth_ids:
            chain_id = str(auth_ids[0])
        elif parsed.get("strand_id"):
            chain_id = str(parsed.get("strand_id") or "").split(",")[0]
        chain = {
            "chain_id": chain_id,
            "auth_asym_id": chain_id,
            "entity_id": parsed.get("entity_id"),
            "sequence": parsed.get("sequence"),
            "length": parsed.get("length"),
            "polymer_type": parsed.get("polymer_type"),
            "description": parsed.get("description"),
            "has_coordinates": False,
        }
        role = classify_chain_role(chain, str(parsed.get("description") or ""))
        chain["role"] = role
        chain["molecule"] = _molecule_from_role(role, str(parsed.get("polymer_type") or ""))
        chains.append(chain)
    if not chains:
        raise ComplexError("No polymer entities to classify.", "PARSING_ERROR")
    guide_map = (
        map_guide_to_rna_chain(guide_sequence, chains)
        if guide_sequence
        else {"status": "UNAVAILABLE", "reason": "No user guide supplied.", "chain_id": ""}
    )
    target_map = (
        map_target_to_dna_chain(target_sequence, chains)
        if target_sequence
        else {"status": "UNAVAILABLE", "reason": "No user target supplied.", "chain_id": ""}
    )
    pam_map = (
        map_pam_on_nontarget(pam, chains)
        if pam or guide_sequence
        else {"status": "UNAVAILABLE", "reason": "No PAM supplied.", "chain_id": ""}
    )
    meta = dict(metadata or {})
    identities = {
        item.get("role"): provenance.sequence_digest(str(item.get("sequence") or ""))
        for item in chains
        if item.get("sequence")
    }
    return provenance.analysis_envelope(
        module="COMPLEX_STRUCTURE",
        payload={
            "structure_id": str(structure_id).upper(),
            "n_chains": len(chains),
            "kind": "experimental",
        },
        status="RETRIEVED",
        algorithm="RCSB Data API polymer_entity (identity only; no Cartesian coordinates)",
        parameters={"structure_id": str(structure_id).upper()},
        source="RCSB PDB",
        accession=str(structure_id).upper(),
    ) | {
        "status": "RETRIEVED",
        "kind": "experimental",
        "coordinates_status": "UNAVAILABLE",
        "structure_id": str(structure_id).upper(),
        "source": "RCSB PDB",
        "method": str(meta.get("method") or ""),
        "resolution_angstrom": meta.get("resolution_angstrom"),
        "chains": chains,
        "atoms": [],
        "n_atoms": 0,
        "synthetic_complex": False,
        "guide_mapping": guide_map,
        "target_mapping": target_map,
        "pam_mapping": pam_map,
        "sequence_identities": identities,
        "mapping_status": {
            "guide": guide_map.get("status"),
            "target": target_map.get("status"),
            "pam": pam_map.get("status"),
        },
        "retrieved_at": provenance.utc_now(),
        "disclaimer": (
            "Chain roles and sequences come from the RCSB polymer_entity API. "
            "Cartesian coordinates are not present in this envelope. HelixScope "
            "did not assemble Cas9, guide and DNA procedurally."
        ),
    }


def complex_provenance(result: Mapping[str, Any]) -> dict:
    """Seis campos de proveniencia de complexo CRISPR.

    Args:
        result: annotate_parsed_complex ou annotate_from_polymer_entities.

    Returns:
        Dict pdb, chains, method, resolution, sequence_identities, mapping_status.

    Raises:
        ComplexError: INVALID_INPUT se PDB ou chains faltarem.
    """
    record = {
        "pdb": str(result.get("structure_id") or ""),
        "chains": [
            {
                "chain_id": item.get("chain_id"),
                "role": item.get("role"),
                "molecule": item.get("molecule"),
            }
            for item in list(result.get("chains") or [])
        ],
        "method": str(result.get("method") or ""),
        "resolution_angstrom": result.get("resolution_angstrom"),
        "sequence_identities": dict(result.get("sequence_identities") or {}),
        "mapping_status": dict(result.get("mapping_status") or {}),
    }
    if not record["pdb"] or not record["chains"]:
        raise ComplexError("Complex provenance requires PDB id and chains.", "INVALID_INPUT")
    return record


def _molecule_from_role(role: str, polymer_type: str) -> str:
    if role in {ROLE_CAS_PROTEIN, ROLE_PROTEIN}:
        return "PROTEIN"
    if role in {ROLE_GUIDE_RNA, ROLE_RNA}:
        return "RNA"
    if role in {ROLE_TARGET_DNA, ROLE_NONTARGET_DNA, ROLE_DNA}:
        return "DNA"
    low = polymer_type.lower()
    if "ribo" in low:
        return "RNA"
    if "deoxy" in low:
        return "DNA"
    if "peptide" in low:
        return "PROTEIN"
    return "UNKNOWN"


def highlight_indices_for_chain_span(
    mappings: Sequence[Mapping[str, Any]],
    chain_id: str,
    start_0based: Optional[int],
    end_0based: Optional[int],
) -> List[int]:
    """Indices de query na cena para um intervalo depositado. Vazio se nao houver mapping.

    Args:
        mappings: Linhas de mapping da cena.
        chain_id: Chain depositada.
        start_0based: Inicio inclusivo no polimero.
        end_0based: Fim exclusivo.

    Returns:
        Lista de query_index_0based. Sem match: lista vazia, nao um highlight inventado.

    Raises:
        Nenhum.
    """
    if start_0based is None or end_0based is None or end_0based <= start_0based:
        return []
    wanted = str(chain_id or "")
    out: List[int] = []
    for item in mappings:
        if str(item.get("chain_id") or "") != wanted:
            continue
        polymer = item.get("chain_index_0based")
        query = item.get("query_index_0based")
        if polymer is None or query is None:
            continue
        if int(start_0based) <= int(polymer) < int(end_0based):
            out.append(int(query))
    return out


def deposited_to_structure_3d_input(
    annotated: Mapping[str, Any],
    *,
    query_sequence: str = "",
    query_molecule: str = "",
) -> dict:
    """Cena 3D a partir de um complexo depositado. Sem XYZ novos.

    Args:
        annotated: Saida de annotate_parsed_complex.
        query_sequence: Sequencia de analise opcional para mapping.
        query_molecule: DNA, RNA ou PROTEIN.

    Returns:
        Contrato structure_3d_input. kind experimental/predicted do ficheiro.

    Raises:
        ComplexError: PARSING_ERROR se nao houver atomos.
    """
    atoms = [dict(item) for item in list(annotated.get("atoms") or [])]
    if not atoms:
        raise ComplexError("Deposited complex has no atoms.", "PARSING_ERROR")
    mappings: List[dict] = []
    scene_index = 0
    query = str(query_sequence or "")
    mol = str(query_molecule or "").strip().upper()
    query_chain = ""
    if query and mol in {"DNA", "RNA", "PROTEIN"}:
        for chain in list(annotated.get("chains") or []):
            polymer = str(chain.get("sequence") or "")
            if polymer and query == polymer:
                query_chain = str(chain.get("chain_id") or "")
                break
    for chain in list(annotated.get("chains") or []):
        cid = str(chain.get("chain_id") or "")
        polymer = str(chain.get("sequence") or "")
        residues = list(chain.get("residues") or [])
        for index, residue in enumerate(residues):
            polymer_index = residue.get("polymer_index_0based")
            if polymer_index is None:
                polymer_index = index
            query_index = None
            query_residue = None
            if query_chain and cid == query_chain and 0 <= int(polymer_index) < len(query):
                query_index = int(polymer_index)
                query_residue = query[query_index]
            elif not query_chain:
                query_index = scene_index
                query_residue = residue.get("one_letter") or (polymer[int(polymer_index)] if int(polymer_index) < len(polymer) else None)
            has = bool(
                (residue.get("ca") or {}).get("x") is not None
                or (residue.get("p") or {}).get("x") is not None
                or (residue.get("c1") or {}).get("x") is not None
            )
            mappings.append(
                {
                    "query_index_0based": query_index,
                    "query_residue": query_residue,
                    "chain_index_0based": polymer_index,
                    "chain_residue": residue.get("one_letter"),
                    "chain_id": cid,
                    "label_seq_id": residue.get("label_seq_id"),
                    "auth_seq_id": residue.get("auth_seq_id"),
                    "insertion_code": residue.get("insertion_code") or "",
                    "comp_id": residue.get("comp_id"),
                    "has_coordinates": bool(residue.get("has_coordinates")) and has,
                    "coordinate_status": "available" if has else "unavailable",
                    "ca": residue.get("ca"),
                    "p": residue.get("p"),
                    "c1": residue.get("c1"),
                    "model": residue.get("model") or chain.get("model") or 1,
                    "molecule_type": chain.get("molecule"),
                    "role": chain.get("role"),
                }
            )
            scene_index += 1
    kind = str(annotated.get("kind") or "experimental")
    n_with = sum(1 for item in mappings if item.get("has_coordinates"))
    return {
        "source": annotated.get("source"),
        "kind": kind,
        "kind_label": annotated.get("kind_label") or protein_structure.structure_kind_label(kind),
        "structure_id": annotated.get("structure_id"),
        "sequence_hash": provenance.sequence_digest(query) if query else "",
        "structure_hash": annotated.get("content_hash"),
        "chain_id": query_chain,
        "model_number": 1,
        "models": sorted({int(item.get("model") or 1) for item in atoms}),
        "chains": list(annotated.get("chains") or []),
        "n_atoms": len(atoms),
        "n_residues": len(mappings),
        "atoms": atoms,
        "residues": [],
        "coordinates": [
            {
                "atom": item.get("atom_name"),
                "element": item.get("element"),
                "residue": item.get("comp_id"),
                "chain": item.get("auth_asym_id") or item.get("label_asym_id"),
                "label_seq_id": item.get("label_seq_id"),
                "auth_seq_id": item.get("auth_seq_id"),
                "insertion_code": item.get("insertion_code") or "",
                "alt_id": item.get("alt_id") or "",
                "x": item.get("x"),
                "y": item.get("y"),
                "z": item.get("z"),
                "occupancy": item.get("occupancy"),
                "b_iso": item.get("b_iso"),
                "model": item.get("model"),
                "group": item.get("group"),
            }
            for item in atoms
        ],
        "mappings": mappings,
        "mapping_status": "mapped" if query_chain else "complex_unmapped_query",
        "mapping_status_note": annotated.get("disclaimer"),
        "coverage_query_with_coordinates": None if not query else (
            sum(1 for item in mappings if item.get("query_index_0based") is not None and item.get("has_coordinates"))
            / max(len(query), 1)
        ),
        "n_with_coordinates": n_with,
        "n_without_coordinates": max(0, len(mappings) - n_with),
        "residues_without_coordinates": [
            item.get("query_index_0based")
            for item in mappings
            if item.get("query_index_0based") is not None and not item.get("has_coordinates")
        ],
        "confidence_data": {},
        "guide_mapping": annotated.get("guide_mapping"),
        "pam_mapping": annotated.get("pam_mapping"),
        "molecule_type": "COMPLEX",
        "metadata": {
            "method": annotated.get("method"),
            "resolution_angstrom": annotated.get("resolution_angstrom"),
            "disclaimer": annotated.get("disclaimer"),
        },
        "provenance": dict(annotated.get("provenance") or {}),
        "renderer": "data-only",
        "status": annotated.get("status") or "RETRIEVED",
        "disclaimer": annotated.get("disclaimer"),
        "synthetic_complex": False,
    }
