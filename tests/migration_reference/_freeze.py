"""Prompt 0 freeze: hashes, env, Streamlit inventory, golden scientific fixtures.

Does not modify scientific algorithms. Does not import Streamlit in the
scientific path. Safe to re-run; overwrites files under tests/migration_reference/.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[2]
REF = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from modules import provenance  # noqa: E402


def _utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def json_safe(value: Any) -> Any:
    """JSON-compatible copy. NumPy/NaN/Path become explicit, not strings-of-science."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, float):
        if math.isnan(value):
            return {"__nonfinite__": "NaN"}
        if math.isinf(value):
            return {"__nonfinite__": "Infinity" if value > 0 else "-Infinity"}
        return float(value)
    if isinstance(value, Path):
        return {"__path__": str(value).replace("\\", "/")}
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return json_safe(value.item())
        if isinstance(value, np.ndarray):
            return {
                "__ndarray__": True,
                "dtype": str(value.dtype),
                "shape": list(value.shape),
                "data": json_safe(value.tolist()),
            }
    except Exception:
        pass
    try:
        import pandas as pd

        if isinstance(value, pd.DataFrame):
            return {
                "__dataframe__": True,
                "columns": list(map(str, value.columns)),
                "records": json_safe(value.to_dict(orient="records")),
            }
    except Exception:
        pass
    if hasattr(value, "__dict__") and type(value).__module__.startswith("Bio"):
        return {"__biopython__": type(value).__name__, "repr": repr(value)[:500]}
    return {"__unserializable__": type(value).__name__}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hash_manifest() -> dict:
    patterns = [
        "modules/*.py",
        "ui/*.py",
        "ui/*.js",
        "app.py",
        "requirements.txt",
        "tests/test_*.py",
        "tests/fixtures/*",
    ]
    files: list[dict] = []
    seen: set[str] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if not path.is_file():
                continue
            rel = path.relative_to(ROOT).as_posix()
            if rel in seen:
                continue
            seen.add(rel)
            files.append(
                {
                    "path": rel,
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                }
            )
    return {
        "generated_at_utc": _utc(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "algorithm": "SHA-256",
        "note": (
            "Evidence of the Prompt 0 freeze. Future extraction may change "
            "bytes; scientific parity is defined by golden fixtures, not this manifest."
        ),
        "n_files": len(files),
        "files": files,
    }


def public_functions(py_path: Path) -> list[str]:
    try:
        tree = ast.parse(py_path.read_text(encoding="utf-8"))
    except SyntaxError:
        return []
    names: list[str] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if not node.name.startswith("_"):
                names.append(node.name)
        elif isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
            names.append(f"class:{node.name}")
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    if not item.name.startswith("_"):
                        names.append(f"{node.name}.{item.name}")
    return names


def classify_module(py_path: Path) -> str:
    text = py_path.read_text(encoding="utf-8", errors="replace")
    has_st = "import streamlit" in text or "from streamlit" in text
    has_sub = "subprocess" in text
    has_http = any(
        token in text
        for token in ("urlopen", "urllib.request", "http.client", "https://")
    )
    has_cache = "@st.cache_data" in text or "@st.cache_resource" in text
    if has_st and any(
        token in text
        for token in ("gc_content", "infer_phylogeny", "pairwise_", "find_guides")
    ):
        return "SCIENCE + STREAMLIT"
    if has_st:
        return "UI ONLY" if "modules" not in py_path.parts else "MIXED / ARCHITECTURAL DEBT"
    if has_sub and has_http:
        return "SCIENCE + SUBPROCESS + REMOTE API"
    if has_sub:
        return "SCIENCE + SUBPROCESS"
    if has_http:
        return "SCIENCE + REMOTE API"
    if has_cache:
        return "SCIENCE + CACHE"
    if py_path.parent.name == "ui" and py_path.name in {
        "charts.py",
        "structure_viewer.py",
        "components.py",
        "tokens.py",
    }:
        return "VISUALIZATION ONLY" if py_path.name in {"charts.py", "structure_viewer.py"} else "UI ONLY"
    return "PURE SCIENCE"


def streamlit_inventory() -> dict:
    files = [ROOT / "app.py", *sorted((ROOT / "ui").glob("*.py")), *sorted((ROOT / "tests").glob("*.py"))]
    usage: dict[str, Any] = {
        "generated_at_utc": _utc(),
        "files": [],
        "session_state_keys": sorted(set(_session_keys())),
        "cache_data_wrappers_in_app_py": [],
        "categories": {
            "UI_ONLY": [],
            "UI_STATE": [],
            "SCIENTIFIC_SESSION_STATE": [],
            "SCIENTIFIC_CACHE": [],
            "FRONTEND_CUSTOM_COMPONENT": [],
            "MIXED / ARCHITECTURAL DEBT": [],
        },
    }
    cache_names: list[str] = []
    for path in files:
        if not path.is_file():
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        imports_streamlit = False
        calls: dict[str, int] = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "streamlit" or alias.name.startswith("streamlit."):
                        imports_streamlit = True
            if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("streamlit"):
                imports_streamlit = True
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "st":
                calls[node.attr] = calls.get(node.attr, 0) + 1
            if isinstance(node, ast.FunctionDef):
                for dec in node.decorator_list:
                    target = dec
                    if isinstance(dec, ast.Call):
                        target = dec.func
                    if (
                        isinstance(target, ast.Attribute)
                        and isinstance(target.value, ast.Name)
                        and target.value.id == "st"
                        and target.attr == "cache_data"
                    ):
                        cache_names.append(node.name)
        if imports_streamlit or calls:
            rel = path.relative_to(ROOT).as_posix()
            usage["files"].append(
                {
                    "path": rel,
                    "imports_streamlit": imports_streamlit,
                    "st_attr_counts": dict(sorted(calls.items(), key=lambda kv: (-kv[1], kv[0]))),
                }
            )
    usage["cache_data_wrappers_in_app_py"] = cache_names
    usage["n_files_importing_streamlit"] = sum(1 for row in usage["files"] if row["imports_streamlit"])
    usage["modules_importing_streamlit"] = []
    for path in sorted((ROOT / "modules").glob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "import streamlit" in text or "from streamlit" in text:
            usage["modules_importing_streamlit"].append(path.name)
    usage["categories"]["SCIENTIFIC_CACHE"] = [
        {"file": "app.py", "function": name, "kind": "st.cache_data wrapper around modules.*"}
        for name in cache_names
    ]
    usage["categories"]["FRONTEND_CUSTOM_COMPONENT"] = [
        {"file": "ui/hand_control.py", "kind": "Streamlit Components v2 + MediaPipe"}
    ]
    usage["categories"]["UI_STATE"] = [
        "helix_active_module",
        "helix_cmd_query",
        "dna_3d_cam_nonce",
        "prot_3d_cam",
    ]
    usage["note"] = (
        "Scientific modules/ do not import Streamlit. Coupling is app.py cache "
        "wrappers and session_state plus ui/ rendering."
    )
    return usage


def _session_keys() -> Iterable[str]:
    text = (ROOT / "app.py").read_text(encoding="utf-8")
    for line in text.splitlines():
        if "session_state[" in line:
            start = line.find("session_state[")
            rest = line[start + len("session_state[") :]
            if rest.startswith(("'", '"')):
                quote = rest[0]
                end = rest.find(quote, 1)
                if end > 1:
                    yield rest[1:end]
    nav = (ROOT / "ui" / "navigation.py").read_text(encoding="utf-8")
    if "ACTIVE_MODULE_KEY" in nav:
        yield "helix_active_module"


def env_snapshot() -> tuple[str, dict]:
    import Bio
    import numpy
    import pandas
    import plotly
    import pytest
    import streamlit

    lines = [
        f"python\t{sys.version.replace(chr(10), ' ')}",
        f"executable\t{sys.executable}",
        f"helixscope_version\t{provenance.HELIXSCOPE_VERSION}",
        f"streamlit\t{streamlit.__version__}",
        f"plotly\t{plotly.__version__}",
        f"biopython\t{Bio.__version__}",
        f"numpy\t{numpy.__version__}",
        f"pandas\t{pandas.__version__}",
        f"pytest\t{pytest.__version__}",
    ]
    try:
        import hypothesis

        lines.append(f"hypothesis\t{hypothesis.__version__}")
    except Exception:
        lines.append("hypothesis\tunavailable")
    try:
        import RNA

        lines.append(f"ViennaRNA_python\t{getattr(RNA, '__version__', 'present')}")
    except Exception:
        lines.append("ViennaRNA_python\tunavailable")
    freeze_path = ROOT / "MIGRATION_PYTHON_ENV.txt"
    freeze_body = "\n".join(lines) + "\n\n# pip freeze (no secrets)\n"
    try:
        import subprocess

        proc = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        freeze_body += proc.stdout or proc.stderr
    except Exception as exc:
        freeze_body += f"# pip freeze failed: {exc}\n"
    tools: dict[str, Any] = {
        "generated_at_utc": _utc(),
        "software_version": provenance.HELIXSCOPE_VERSION,
        "python": sys.version.split()[0],
        "packages": {
            "streamlit": streamlit.__version__,
            "plotly": plotly.__version__,
            "biopython": Bio.__version__,
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
        },
        "env_var_names_no_values": sorted(
            {
                "NCBI_API_KEY",
                "ENTREZ_API_KEY",
                "HELIXSCOPE_TOOLS_DIR",
                "HELIXSCOPE_IQTREE",
                "HELIXSCOPE_FASTTREE",
                "HELIXSCOPE_DSSP",
                "HELIXSCOPE_STRIDE",
                "HELIXSCOPE_BLASTN",
                "HELIXSCOPE_BLAST_DB",
                "HELIXSCOPE_MAKEBLASTDB",
                "HELIXSCOPE_RNAFOLD",
                "HELIXSCOPE_CAS_OFFINDER",
                "HELIXSCOPE_REFERENCE_DIR",
                "HELIXSCOPE_JOBS_DIR",
                "HELIXSCOPE_STRUCTURE_CACHE",
                "HELIXSCOPE_CONTACT_EMAIL",
                "HELIXSCOPE_LIVE_GENOME",
                "HELIXSCOPE_OPENCL_PROBE_WMI",
            }
        ),
        "env_present_names_only": sorted(
            name
            for name in (
                "NCBI_API_KEY",
                "ENTREZ_API_KEY",
                "HELIXSCOPE_TOOLS_DIR",
                "HELIXSCOPE_IQTREE",
                "HELIXSCOPE_FASTTREE",
                "HELIXSCOPE_DSSP",
                "HELIXSCOPE_BLAST_DB",
                "HELIXSCOPE_CAS_OFFINDER",
                "HELIXSCOPE_REFERENCE_DIR",
                "HELIXSCOPE_JOBS_DIR",
                "HELIXSCOPE_STRUCTURE_CACHE",
                "HELIXSCOPE_CONTACT_EMAIL",
                "HELIXSCOPE_LIVE_GENOME",
            )
            if (os.environ.get(name) or "").strip()
        ),
    }
    freeze_path.write_text(freeze_body, encoding="utf-8")
    return freeze_body, tools


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def golden_dna() -> dict:
    from modules import dna_analysis, provenance as prov

    cases = {
        "short_valid": "ATGCATGCATGC",
        "mixed_composition": "ACGTACGTAAAATTTTGGGGCCCC",
        "high_gc": "GCGCGCGCGCGCGCGC",
        "low_gc": "ATATATATATATATAT",
        "ambiguous_n": "ATGCNNNNATGC",
        "all_n": "NNNN",
        "scale_50k": ("ATGC" * 12_500),
    }
    out: dict[str, Any] = {
        "parity_class_notes": {
            "sequence": "EXACT",
            "gc_content": "FLOAT_TOLERANCE abs=1e-12 relative=0 (Python float)",
            "integer_counts": "EXACT",
            "NaN": "must remain non-finite / Unavailable, never 0",
        },
        "cases": {},
    }
    for name, seq in cases.items():
        info = dna_analysis.validate_sequence(seq)
        cleaned = str(info.get("sequence") or "")
        gc = dna_analysis.gc_content(cleaned) if cleaned else float("nan")
        record: dict[str, Any] = {
            "input": seq if name != "scale_50k" else {"repeat": "ATGC", "n_repeats": 12500, "length": 50000},
            "validation": {
                "is_valid": info.get("is_valid"),
                "type": info.get("type"),
                "length": info.get("length") or len(cleaned),
                "rejection_reason": info.get("rejection_reason"),
            },
            "sequence_hash": prov.sequence_digest(cleaned) if cleaned else "",
            "gc_content": gc,
            "at_content": dna_analysis.at_content(cleaned) if cleaned else float("nan"),
            "composition": dna_analysis.nucleotide_composition(cleaned) if cleaned else {},
            "reverse_complement": (
                dna_analysis.reverse_complement(cleaned) if cleaned else None
            ),
            "shannon_entropy": dna_analysis.shannon_entropy(cleaned) if cleaned else float("nan"),
        }
        if name != "scale_50k" and cleaned and set(cleaned) <= set("ACGT"):
            record["melting_temperature_c"] = dna_analysis.melting_temperature(cleaned)
            record["molecular_weight_dna"] = dna_analysis.molecular_weight(cleaned, "DNA")
            record["tm_report"] = dna_analysis.melting_temperature_report(cleaned)
        if name == "short_valid":
            record["orfs_min_12"] = dna_analysis.find_orfs(cleaned, min_length=12)
        out["cases"][name] = record
    return out


def golden_rna() -> dict:
    from modules import dna_analysis, rna_analysis, rna_folding

    dna = "ATGCATGCATGC"
    rna = rna_analysis.transcribe(dna)
    coding = "AUG" + "GCU" * 20 + "UAA"
    fold_seq = "GCGCGCAAAAGCGCGC"
    let7 = ""
    let7_path = ROOT / "tests" / "fixtures" / "rna_hsa_let7a.fa"
    if let7_path.is_file():
        parsed = dna_analysis.parse_sequence_payload(let7_path.read_text(encoding="utf-8"))
        let7 = str(parsed.get("sequence") or "")
    fold: dict[str, Any] = {}
    try:
        fold = rna_folding.fold_rna(fold_seq, identifier="prompt0_hairpin", backend="viennarna_python")
        for drop in ("computed_timestamp", "retrieval_timestamp"):
            fold.pop(drop, None)
            if isinstance(fold.get("provenance"), dict):
                fold["provenance"].pop(drop, None)
    except Exception as exc:
        fold = {"status": "UNAVAILABLE", "error_type": type(exc).__name__, "message": str(exc)}
    cai = None
    enc = None
    try:
        cai = rna_analysis.codon_adaptation_index_report(coding, "human")
        enc = rna_analysis.effective_number_of_codons_report(coding)
    except Exception as exc:
        cai = {"status": "ERROR", "message": str(exc)}
    return {
        "transcription": {"dna": dna, "rna": rna},
        "composition": {
            "rna_sequence": rna,
            "shannon_entropy": rna_analysis.shannon_entropy(rna),
        },
        "codon_fixture": {
            "coding_rna": coding,
            "cai_report": cai,
            "enc_report": enc,
        },
        "fold_python_viennarna": fold,
        "let7_length": len(let7),
        "let7_not_folded_here": True,
        "note": "ViennaRNA Python MFE is PREDICTED. RNAfold CLI is a different backend.",
    }


def golden_protein() -> dict:
    from modules import protein_analysis

    # Crambin-like short standard protein from 1CRN sequence context; use BLAST tiny protein.
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    report = protein_analysis.physicochemical_report(seq, ph=7.0)
    return {
        "sequence": seq,
        "length": len(seq),
        "composition": protein_analysis.amino_acid_composition(seq),
        "physicochemical_report": report,
        "ss_prediction_status": protein_analysis.secondary_structure_availability(),
        "note": "Sequence SS predictors remain UNAVAILABLE. DSSP is coordinate assignment, not this report.",
    }


def golden_alignment() -> dict:
    from modules import alignment

    a = "ACGTACGTACGT"
    b = "ACGTACGGACGT"
    return {
        "global": alignment.pairwise_global(a, b),
        "local": alignment.pairwise_local(a, b),
        "parameters": {
            "match": alignment.MATCH_SCORE,
            "mismatch": alignment.MISMATCH_SCORE,
            "open_gap": alignment.OPEN_GAP_SCORE,
            "extend_gap": alignment.EXTEND_GAP_SCORE,
        },
    }


def golden_motif() -> dict:
    from modules import motif_search

    seq = "ACGTACGTACGT"
    return {
        "exact": motif_search.find_motif(seq, "ACGT"),
        "iupac": motif_search.find_motif(seq, "NNGN"),
        "coordinate_convention": "start 0-based inclusive, end exclusive, on the informed strand",
    }


def golden_phylogeny() -> dict:
    from modules import msa, phylogeny

    identifiers = ["seq_a", "seq_b", "seq_c", "seq_d"]
    aligned = [
        "ACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTTTTT",
        "TTTTACGTACGTACGTACGT",
        "ACGTACGTAAAAACGTACGT",
    ]
    members = [
        msa.build_collection_member(
            sequence=row.replace("-", ""),
            identifier=ident,
            source="user input",
            molecule="DNA",
        )
        for ident, row in zip(identifiers, aligned)
    ]
    raw = "".join(f">{i}\n{row}\n" for i, row in zip(identifiers, aligned))
    msa_result = msa.build_msa_result(
        raw_alignment=raw,
        members=members,
        tool="fixture",
        tool_version="prompt0",
        method="prealigned fixture",
        parameters={"source": "tests/phase13_streamlit_flow.PREALIGNED"},
    )
    payload: dict[str, Any] = {
        "msa_alignment_hash": msa_result.get("alignment_hash"),
        "n_sequences": msa_result.get("n_sequences"),
        "n_columns": msa_result.get("n_columns"),
        "trees": {},
    }
    methods = [
        ("nj", phylogeny.METHOD_NJ, {"distance_model": phylogeny.DISTANCE_P}),
        ("upgma", phylogeny.METHOD_UPGMA, {"distance_model": phylogeny.DISTANCE_P}),
        ("iqtree", phylogeny.METHOD_IQTREE, {}),
        ("fasttree", phylogeny.METHOD_FASTTREE, {}),
    ]
    for key, method, extra in methods:
        try:
            tree = phylogeny.infer_phylogeny(msa_result, method=method, **extra)
            slim = {
                "status": tree.get("status"),
                "method": tree.get("method"),
                "tool": tree.get("tool"),
                "tool_version": tree.get("tool_version"),
                "n_leaves": tree.get("n_leaves"),
                "newick": tree.get("newick"),
                "model": tree.get("model") or tree.get("selected_model") or extra.get("distance_model"),
                "support": tree.get("support"),
                "not_taxonomic_tree": tree.get("not_taxonomic_tree"),
                "alignment_hash": tree.get("alignment_hash"),
                "leaves": [
                    {"tree_id": leaf.get("tree_id"), "sequence_hash": leaf.get("sequence_hash")}
                    for leaf in (tree.get("leaves") or [])
                ],
                "disclaimer": tree.get("disclaimer"),
            }
            payload["trees"][key] = slim
        except Exception as exc:
            payload["trees"][key] = {
                "status": "ERROR",
                "error_type": type(exc).__name__,
                "message": str(exc),
            }
    payload["semantic_parity"] = (
        "Newick sibling order may differ across engines. Compare leaf sets, "
        "status, method, and engine-specific model/support semantics. Do not "
        "require NJ == IQ-TREE topology."
    )
    return payload


def golden_crispr() -> dict:
    from modules import crispr, crispr_reference

    guide = "ACGTACGTACGTACGTACGT"
    target = guide + "AGG"
    found = crispr.find_guides(target, "SpCas9")
    availability = crispr.model_availability()
    ref_text = (ROOT / "tests" / "fixtures" / "helixscope_test_reference.fa").read_text(encoding="utf-8")
    reference = crispr_reference.load_reference_from_text(
        ref_text,
        source="TEST REFERENCE fixture",
        organism="synthetic",
        assembly="HELIXSCOPE_TEST_REF",
        version="1",
        accession="HELIXSCOPE_TEST_REF",
        upload_filename="helixscope_test_reference.fa",
    )
    return {
        "scope": "TEST REFERENCE / pasted sequence. Not genome-wide. Not GRCh38.",
        "guide_scan": found,
        "model_availability": availability,
        "reference_identity": {
            "assembly": reference.get("assembly") or reference.get("assembly_id"),
            "status": reference.get("status"),
            "n_contigs": reference.get("n_contigs") or len(reference.get("records") or []),
        },
        "on_target_doench_ruleset2_available": bool(
            (availability.get("on_target_doench_ruleset2") or {}).get("available")
        ),
        "on_target_deephf_available": bool(
            (availability.get("on_target_deephf") or {}).get("available")
        ),
    }


def golden_variant() -> dict:
    from modules import variant_core

    v38 = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh38.p14")
    v37 = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh37")
    return {
        "note": "Identity only. No VEP/ClinVar live call in this golden. BRCA1 locus example is coordinate identity.",
        "grch38": {
            "identity_hash": v38.get("identity_hash"),
            "contig": v38.get("contig"),
            "position_0based": v38.get("position_0based"),
            "position_1based": v38.get("position_1based"),
            "ref": v38.get("ref"),
            "alt": v38.get("alt"),
            "assembly": v38.get("assembly"),
            "status": v38.get("status"),
            "effect_status": v38.get("effect_status"),
        },
        "same_text_different_assembly_identity_differs": v38.get("identity_hash") != v37.get("identity_hash"),
    }


def golden_structures() -> dict:
    fixtures = ROOT / "tests" / "fixtures"
    records = {}
    for name in ("1BNA.cif", "1RNA.cif", "1CRN.cif"):
        path = fixtures / name
        if path.is_file():
            records[name] = {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "source": "tests/fixtures (vendored)",
                "kind": "experimental fixture",
            }
    cache = ROOT / "data" / "structure_cache"
    four = None
    if cache.is_dir():
        for path in cache.glob("*"):
            if path.suffix.lower() in {".cif", ".pdb"} and path.stat().st_size > 1_000_000:
                four = {
                    "path_basename": path.name,
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                    "note": "Local cache of RCSB download if present. Mapping of 4UN3 to CRISPR guide remains UNCERTAIN.",
                }
                break
    return {
        "fixtures": records,
        "4UN3_cache": four,
        "mapping_law": "UNCERTAIN mapping must not be promoted to EXACT by a new frontend.",
        "kinds": ["experimental", "predicted", "illustrative", "unavailable"],
    }


def golden_compare() -> dict:
    live = ROOT / "data" / "jobs" / "engine_live.json"
    payload: dict[str, Any] = {
        "pair": "8HSK vs 8HSF",
        "method": "fatcat-rigid / jFATCAT-rigid",
        "source": "data/jobs/engine_live.json REMOTE_VALIDATED snapshot, not a new RCSB call",
        "expected_semantics": {
            "n_aligned_ca_pairs": 20,
            "block_rmsd_angstrom": 0.8,
            "global_rmsd_angstrom": 0.99,
            "tm_score": 0.23,
            "coverage_note": "95% / 100% where the current path supplies it",
        },
        "parity": "FLOAT_TOLERANCE RMSD abs=0.05 A (RCSB reported precision 0.01-0.1). Pair count EXACT.",
        "live_record": None,
    }
    if live.is_file():
        data = json.loads(live.read_text(encoding="utf-8"))
        payload["live_record"] = data.get("RCSB Alignment API")
    fixture = ROOT / "tests" / "fixtures" / "rcsb_alignment_complete.json"
    if fixture.is_file():
        from modules import structure_alignment

        parsed = structure_alignment.parse_alignment_payload(
            json.loads(fixture.read_text(encoding="utf-8"))
        )
        payload["fixture_parser"] = {
            "n_pairs": parsed.get("n_aligned_residue_pairs") or parsed.get("n_pairs"),
            "status": parsed.get("status"),
            "method": parsed.get("method"),
        }
    return payload


def golden_evidence() -> dict:
    from modules import evidence_workspace

    items = [
        evidence_workspace.evidence_item(
            field="Variant A most_severe_consequence",
            value="missense_variant",
            source="Ensembl VEP",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
        evidence_workspace.evidence_item(
            field="Variant B most_severe_consequence",
            value="synonymous_variant",
            source="NCBI ClinVar",
            evidence_status="RETRIEVED",
            retrieved_at_utc="2026-08-29T00:00:00Z",
        ),
    ]
    conflicts = evidence_workspace.source_conflicts(items)
    pack = evidence_workspace.export_pack(
        kind="variants",
        inputs={"example": "conflict"},
        items=items,
    )
    return {
        "items": items,
        "conflicts": conflicts,
        "confidence_score": pack.get("confidence_score"),
        "policy": "Conflicts are listed. Majority vote does not resolve them. confidence_score is None.",
    }


def golden_blast() -> dict:
    from modules import blast_search

    avail = blast_search.local_blast_availability()
    result: dict[str, Any] = {
        "availability": {
            "executables_detected": bool(avail.get("executables_detected")),
            "version": avail.get("version"),
            "not_nt_nr": True,
        },
        "run": None,
    }
    if not avail.get("executables_detected") or not blast_search.detect_makeblastdb().get("available"):
        result["run"] = {"status": "NOT_INSTALLED", "note": "BLAST+ not used for this freeze run."}
        return result
    fasta = (ROOT / "tests" / "fixtures" / "blast_tiny_nucl.fa").read_text(encoding="utf-8")
    query = "ATGCTAGTCGGATCCTGAATGCGTACGACTAG"
    original = os.environ.get("HELIXSCOPE_BLAST_DB")
    try:
        with tempfile.TemporaryDirectory() as tmp:
            prefix = str(Path(tmp) / "tiny_nucl")
            built = blast_search.build_declared_blast_database(
                fasta_text=fasta,
                dbtype="nucl",
                prefix=prefix,
                title="helixscope_tiny_nucl",
            )
            os.environ["HELIXSCOPE_BLAST_DB"] = prefix
            hit = blast_search.run_local_blast(program="blastn", query=query, molecule="DNA")
            result["run"] = {
                "db_not_nt": built.get("not_nt"),
                "backend": hit.get("backend"),
                "status": hit.get("status"),
                "n_hits": hit.get("n_hits"),
                "not_remote_ncbi": hit.get("not_remote_ncbi"),
            }
    except Exception as exc:
        result["run"] = {"status": "ERROR", "error_type": type(exc).__name__, "message": str(exc)}
    finally:
        if original is None:
            os.environ.pop("HELIXSCOPE_BLAST_DB", None)
        else:
            os.environ["HELIXSCOPE_BLAST_DB"] = original
    return result


def function_inventory() -> dict:
    modules_dir = ROOT / "modules"
    rows = []
    n_public = 0
    for path in sorted(modules_dir.glob("*.py")):
        if path.name == "__init__.py":
            continue
        fns = public_functions(path)
        n_public += len([n for n in fns if not n.startswith("class:")])
        rows.append(
            {
                "path": path.relative_to(ROOT).as_posix(),
                "classification": classify_module(path),
                "n_public_symbols": len(fns),
                "public": fns,
            }
        )
    ui_rows = []
    for path in sorted((ROOT / "ui").glob("*.py")):
        ui_rows.append(
            {
                "path": path.relative_to(ROOT).as_posix(),
                "classification": classify_module(path),
                "imports_streamlit": "import streamlit" in path.read_text(encoding="utf-8"),
            }
        )
    return {
        "n_module_files": len(rows),
        "n_public_module_symbols": n_public,
        "n_module_files_importing_streamlit": 0,
        "modules": rows,
        "ui": ui_rows,
        "app_py": "SCIENCE + STREAMLIT cache wrappers + UI composition",
    }


def detect_tools() -> dict:
    from modules import tool_registry

    snap = tool_registry.collect_tool_snapshot(refresh=True)
    slim = {}
    for name, row in (snap.get("tools") or {}).items():
        slim[name] = {
            "status": row.get("status"),
            "version": row.get("version"),
            "available": row.get("available"),
            "source": row.get("source"),
            "engine_location": row.get("engine_location"),
            "live": row.get("live"),
        }
    return slim


def main() -> int:
    REF.mkdir(parents=True, exist_ok=True)
    hashes = hash_manifest()
    write_json(ROOT / "MIGRATION_REFERENCE_HASHES.json", hashes)
    _, tools_pkg = env_snapshot()
    tools_pkg["engines"] = detect_tools()
    write_json(ROOT / "MIGRATION_TOOL_VERSIONS.json", tools_pkg)
    write_json(ROOT / "MIGRATION_STREAMLIT_DEPENDENCIES.json", streamlit_inventory())
    write_json(REF / "function_inventory.json", function_inventory())
    write_json(REF / "dna.json", golden_dna())
    write_json(REF / "rna.json", golden_rna())
    write_json(REF / "protein.json", golden_protein())
    write_json(REF / "alignment.json", golden_alignment())
    write_json(REF / "motif.json", golden_motif())
    write_json(REF / "phylogeny.json", golden_phylogeny())
    write_json(REF / "crispr.json", golden_crispr())
    write_json(REF / "variant.json", golden_variant())
    write_json(REF / "structures.json", golden_structures())
    write_json(REF / "compare.json", golden_compare())
    write_json(REF / "evidence.json", golden_evidence())
    write_json(REF / "blast.json", golden_blast())
    write_json(
        REF / "parity_classes.json",
        {
            "EXACT": [
                "sequence",
                "identifier",
                "status",
                "integer count",
                "coordinate integers",
                "hash",
                "accession",
                "boolean semantic flags",
                "Newick leaf set (not sibling order)",
            ],
            "FLOAT_TOLERANCE": {
                "gc_fraction": {"abs": 1e-12, "justification": "IEEE-754 identity of the same Python formula"},
                "Tm_C": {"abs": 1e-9, "justification": "same SantaLucia/Wallace implementation"},
                "MW": {"abs": 1e-6, "unit": "Da or kDa as stored", "justification": "ProtParam float"},
                "RMSD_A": {"abs": 0.05, "justification": "RCSB/US-align reported precision; do not use ±10%"},
                "SASA_A2": {"abs": 0.5, "justification": "Shrake-Rupley probe sampling"},
                "MFE_kcal_mol": {"abs": 1e-4, "justification": "ViennaRNA Python float32/float64 mix observed on 2.7.2"},
                "identity_pct": {"abs": 0.01, "justification": "stored at 2 decimal places"},
            },
            "STRUCTURAL_ARRAY_PARITY": "mmCIF coordinates identical when the same file hash is used",
            "SEMANTIC_PARITY": ["Newick sibling order", "JSON key order", "dict insertion order"],
        },
    )
    print("Prompt 0 freeze written", REF)
    print("files", hashes["n_files"], "version", provenance.HELIXSCOPE_VERSION)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
