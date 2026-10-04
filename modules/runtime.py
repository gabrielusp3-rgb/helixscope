"""Contratos de API do processo Streamlit.

No Windows o Streamlit reexecuta app.py no mesmo interpretador. Imports
`from modules import x` devolvem o objeto antigo em sys.modules. Este modulo
recarrega UM modulo apenas quando um atributo exigido falta, no maximo uma
vez por (modulo, versao esperada). Se o atributo continuar ausente, o
resultado e needs_restart — nao um reload em massa a cada rerun.

Nenhuma funcao importa Streamlit.
"""

from __future__ import annotations

import importlib
import sys
from typing import Any, Dict, Iterable, Optional, Tuple

REQUIRED_APIS: Dict[str, Tuple[str, ...]] = {
    "ui.structure_viewer": (
        "parse_plotly_selection",
        "parse_plotly_selection_indices",
        "figure_from_scene",
        "figure_from_superposition",
        "STRUCTURE_VIEWER_API",
    ),
    "modules.provenance": (
        "HELIXSCOPE_VERSION",
        "hashes_match",
        "sequence_digest",
        "analysis_envelope",
    ),
    "modules.thermodynamics": ("santalucia_nearest_neighbor_tm",),
    "modules.tool_registry": ("collect_tool_snapshot",),
    "modules.phylogeny": ("infer_phylogeny", "detect_iqtree", "METHOD_IQTREE"),
    "modules.protein_structure": (
        "assign_secondary_structure_dssp",
        "shrake_rupley_sasa",
        "map_dssp_assignments",
        "interchain_contacts",
    ),
    "modules.msa": ("build_msa_result", "import_prealigned_fasta"),
    "modules.blast_search": (
        "local_blast_availability",
        "run_local_blast",
        "ncbi_blast_availability",
    ),
    "modules.genome_fasta": (
        "build_fai_for_fasta",
        "fetch_sequence",
        "fai_is_current",
    ),
    "modules.genome_store": (
        "ready_record",
        "list_local_references",
        "install_test_reference",
    ),
    "modules.genome_jobs": (
        "submit_casoffinder_job",
        "load_status",
        "load_result",
        "request_cancel",
        "full_genome_search_preflight",
    ),
    "modules.resource_admission": ("admit_reference_search", "system_memory"),
    "modules.opencl_runtime": ("diagnose",),
    "modules.genome_annotation": ("parse_gff_text", "classify_hit_context"),
    "modules.rcsb_alignment": (
        "align_pairwise",
        "supported_methods",
        "url_is_allowed",
        "network_disclosure",
    ),
    "modules.structure_alignment": ("parse_alignment_payload", "lookup_correspondence"),
    "modules.structure_superposition": (
        "apply_column_major_4x4",
        "independent_aligned_ca_rmsd",
        "superposition_bundle",
    ),
    "modules.comparative": (
        "compare_structures_remote",
        "compare_structures_usalign",
        "compare_variants",
        "compare_guides",
    ),
    "modules.evidence_workspace": ("evidence_item", "source_conflicts", "export_pack"),
    "modules.usalign": (
        "detect_usalign",
        "align_structure_files",
        "parse_usalign_stdout",
    ),
    "modules.explain": ("explain", "display_number"),
}

_ENSURED_FOR_VERSION: str = ""
_LAST_REPORT: Dict[str, Any] = {}


def load_module(module_name: str) -> Any:
    """Importa um modulo pelo nome pontilhado.

    Args:
        module_name: Nome em sys.modules (ex. ui.structure_viewer).

    Returns:
        Objeto modulo.

    Raises:
        ImportError: Se o modulo nao existir.
    """
    return importlib.import_module(module_name)


def missing_attributes(module_obj: Any, names: Iterable[str]) -> list[str]:
    """Lista atributos publicos ausentes.

    Args:
        module_obj: Modulo ja importado.
        names: Nomes exigidos.

    Returns:
        Nomes que getattr devolve None e que nao existem.

    Raises:
        Nenhum.
    """
    missing: list[str] = []
    for name in names:
        if not hasattr(module_obj, name):
            missing.append(str(name))
    return missing


def reload_if_missing(module_name: str, names: Iterable[str]) -> dict:
    """Recarrega o modulo uma vez se faltar algum atributo.

    Args:
        module_name: Nome pontilhado.
        names: Atributos exigidos.

    Returns:
        Dict module, missing_before, missing_after, reloaded.

    Raises:
        ImportError: Se o import falhar.
    """
    module_obj = load_module(module_name)
    before = missing_attributes(module_obj, names)
    reloaded = False
    if before:
        module_obj = importlib.reload(module_obj)
        reloaded = True
    after = missing_attributes(module_obj, names)
    return {
        "module": module_name,
        "object": module_obj,
        "missing_before": before,
        "missing_after": after,
        "reloaded": reloaded,
    }


def ensure_runtime(*, expected_version: str) -> dict:
    """Garante contratos de API. Nao recarrega dezenas de modulos por rerun.

    Args:
        expected_version: Versao do codigo em app.py / provenance.

    Returns:
        Dict ok, needs_restart, expected_version, loaded_version, reloaded,
        missing, restart_reason.

    Raises:
        Nenhum. Falhas de import entram em missing.
    """
    global _ENSURED_FOR_VERSION, _LAST_REPORT
    version = str(expected_version or "").strip()
    if _ENSURED_FOR_VERSION == version and _LAST_REPORT.get("ok"):
        return dict(_LAST_REPORT)

    reloaded: list[str] = []
    missing: list[str] = []
    loaded_version = ""
    for module_name, attrs in REQUIRED_APIS.items():
        try:
            result = reload_if_missing(module_name, attrs)
        except ImportError as exc:
            missing.append(f"{module_name}: ImportError {exc}")
            continue
        if result["reloaded"]:
            reloaded.append(module_name)
        for attr in result["missing_after"]:
            missing.append(f"{module_name}.{attr}")
        if module_name == "modules.provenance":
            loaded_version = str(getattr(result["object"], "HELIXSCOPE_VERSION", "") or "")

    version_mismatch = bool(loaded_version) and bool(version) and loaded_version != version
    needs_restart = bool(missing) or version_mismatch
    reason = ""
    if missing:
        reason = (
            "HelixScope APIs are missing in this process (stale Python module). "
            "Stop the app and run: python -m streamlit run app.py"
        )
    elif version_mismatch:
        reason = (
            f"Loaded provenance version {loaded_version} differs from "
            f"expected {version}. Stop Streamlit and restart the process."
        )
    report = {
        "ok": not needs_restart,
        "needs_restart": needs_restart,
        "expected_version": version,
        "loaded_version": loaded_version,
        "reloaded": reloaded,
        "missing": missing,
        "restart_reason": reason,
        "scan": "attribute contracts only; no mass importlib.reload on every rerun",
    }
    if report["ok"]:
        _ENSURED_FOR_VERSION = version
    _LAST_REPORT = dict(report)
    return dict(report)


def require_attr(module_obj: Any, name: str) -> Any:
    """Devolve um atributo, recarregando o modulo uma vez se estiver ausente.

    Args:
        module_obj: Modulo.
        name: Atributo publico.

    Returns:
        O atributo, ou None se continuar ausente apos um reload.

    Raises:
        Nenhum.
    """
    if module_obj is None:
        return None
    attr = getattr(module_obj, name, None)
    if attr is not None:
        return attr
    reloaded = importlib.reload(module_obj)
    return getattr(reloaded, name, None)


def reset_runtime_cache() -> None:
    """Limpa o atalho in-process (testes).

    Args:
        Nenhum.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    global _ENSURED_FOR_VERSION, _LAST_REPORT
    _ENSURED_FOR_VERSION = ""
    _LAST_REPORT = {}


def module_is_current(module_name: str, attr: str) -> bool:
    """True se sys.modules tem o modulo com o atributo.

    Args:
        module_name: Nome pontilhado.
        attr: Atributo exigido.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    module_obj = sys.modules.get(module_name)
    if module_obj is None:
        return False
    return hasattr(module_obj, attr)
