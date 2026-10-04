"""Testes de regressao das correcoes de seguranca e limites tecnicos.

Nao disparam consultas de rede. Nao tentam explorar servicos de terceiros.
"""

import pytest

from modules import alignment, blast_search, crispr, crispr_casoffinder, dna_analysis, motif_search, msa, ncbi_fetch, phylogeny, protein_structure, rna_folding
from ui import components


def test_html_components_escape_untrusted_text():
    payload = '<img src=x onerror="alert(1)">'
    card = components.glass_card(payload, "safe-body")
    chip = components.metric_chip(payload, payload, payload)
    pill = components.badge(payload, "warn")
    viewer = components.sequence_display(payload + "ACGT", "DNA")
    from ui.workspace import metric_grid_html

    grid = metric_grid_html([(payload, payload, payload)])
    for html in (card, chip, pill, viewer, grid):
        assert "<IMG" not in html
        assert "<img" not in html
        assert "&lt;" in html
    status = components.status_badge("COMPUTED")
    assert "COMPUTED" in status
    escaped = components.status_badge('<img src=x onerror="alert(1)">')
    assert "<img" not in escaped
    assert "&lt;" in escaped


def test_sequence_display_truncates_viewer_only():
    seq = "A" * (components.MAX_SEQUENCE_DISPLAY_RESIDUES + 50)
    html = components.sequence_display(seq, "DNA")
    assert "Showing residues" in html
    assert f"{components.MAX_SEQUENCE_DISPLAY_RESIDUES + 50:,}" in html


def test_safe_download_filename_blocks_path_traversal():
    name = components.safe_download_filename(r"..\..\etc\passwd")
    assert ".." not in name
    assert "\\" not in name
    assert "/" not in name
    assert name == "etc_passwd"
    assert components.safe_download_filename("<script>alert(1)</script>") == "script_alert_1_script"
    assert components.safe_download_filename("...") == "download"


def test_validate_rejects_oversized_raw_and_residue_input():
    raw = "A" * (dna_analysis.MAX_RAW_INPUT_CHARS + 1)
    raw_result = dna_analysis.validate_sequence(raw)
    assert raw_result["is_valid"] is False
    assert "technical limit" in raw_result["rejection_reason"]

    residues = "A" * (dna_analysis.MAX_INPUT_RESIDUES + 1)
    mol = dna_analysis.validate_for_molecule(residues, "DNA")
    assert mol["is_valid"] is False
    assert "technical limit" in mol["rejection_reason"]
    short = dna_analysis.validate_for_molecule("ATGCATGC", "DNA")
    assert short["is_valid"] is True


def test_alignment_and_dotplot_reject_oversized_pairs():
    long_seq = "A" * (alignment.MAX_ALIGN_LENGTH + 1)
    with pytest.raises(ValueError, match="limited"):
        alignment.pairwise_global(long_seq, "ACGT")
    with pytest.raises(ValueError, match="limited"):
        alignment.pairwise_local(long_seq, "ACGT")
    huge_dot = "A" * (alignment.MAX_DOTPLOT_LENGTH + 1)
    with pytest.raises(ValueError, match="limited"):
        alignment.dotplot_matrix(huge_dot, "ACGT")


def test_ncbi_rejects_oversized_identifier_and_invalid_email():
    huge = "A" * (ncbi_fetch.MAX_IDENTIFIER_LENGTH + 1)
    info = ncbi_fetch.classify_identifier(huge)
    assert info["kind"] == "invalid"
    gene = ncbi_fetch.classify_identifier("1" * (ncbi_fetch.MAX_GENE_ID_DIGITS + 1))
    assert gene["kind"] == "invalid"
    with pytest.raises(ValueError, match="e-mail"):
        ncbi_fetch.fetch_by_accession("NM_000518.5", "not-an-email")
    with pytest.raises(ValueError, match="e-mail"):
        ncbi_fetch.fetch_by_accession("NM_000518.5", "user@")
    assert ncbi_fetch._valid_entrez_email("user@example.com") is True


def test_entrez_urlopen_sets_timeout(monkeypatch):
    captured = {}

    def fake_urlopen(url, data=None, timeout=None, *, context=None):
        captured["timeout"] = timeout
        class Handle:
            pass

        return Handle()

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    ncbi_fetch._entrez_urlopen("https://example.invalid/entrez")
    assert captured["timeout"] == ncbi_fetch.ENTREZ_TIMEOUT_S


def test_crispr_off_target_scan_rejects_oversized_target():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    huge = "A" * (crispr.MAX_OFF_TARGET_SCAN_NT + 1)
    with pytest.raises(ValueError, match="off-target"):
        crispr.find_off_targets(guide, huge, 3)
    found = [{"guide_sequence": guide, "pam_sequence": "AGG", "position": 1, "strand": "+"}]
    with pytest.raises(ValueError, match="off-target"):
        crispr.evaluate_guides(found, huge, "SpCas9", run_off_target=True)


def test_crispr_report_index_is_not_a_spreadsheet_formula():
    report = crispr.generate_report(
        [
            {
                "rank": 1,
                "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
                "pam_sequence": "AGG",
                "position": 1,
                "strand": "+",
                "gc_content": 50.0,
                "doench_score": 0.7,
            }
        ],
        "=cmd|'/c calc'!A0",
    )
    assert not str(report.index.name).startswith("=")


def test_motif_pattern_rejects_extreme_length():
    with pytest.raises(ValueError, match="excede"):
        motif_search.find_motif("ACGT", "N" * (motif_search.MAX_MOTIF_PATTERN_LENGTH + 1))


def test_kmer_and_alignment_caps_fail_closed():
    with pytest.raises(ValueError):
        dna_analysis.kmer_counts("ACGT", k=1_000_000)
    with pytest.raises(ValueError, match="limited"):
        alignment.pairwise_global("A" * (alignment.MAX_ALIGN_LENGTH + 1), "ACGT")
    with pytest.raises(ValueError, match="limited"):
        alignment.dotplot_matrix("A" * (alignment.MAX_DOTPLOT_LENGTH + 1), "ACGT")


def test_crispr_pam_scan_rejects_oversized_sequence():
    huge = "A" * (crispr.MAX_GUIDE_SCAN_NT + 1)
    with pytest.raises(ValueError, match="limited"):
        crispr.find_guides(huge, "SpCas9")


def test_blast_rid_and_endpoint_are_not_user_urls():
    assert blast_search.BLAST_ENDPOINT == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    with pytest.raises(blast_search.BlastError) as rid:
        blast_search.poll_search(
            {
                "rid": "http://127.0.0.1/Blast.cgi",
                "status": "WAITING",
                "submitted_monotonic": 0.0,
            },
            "user@example.com",
            now=10_000.0,
        )
    assert rid.value.category == "INVALID_INPUT"
    with pytest.raises(blast_search.BlastError) as xml:
        blast_search.parse_blast_xml("A" * (blast_search.MAX_XML_BYTES + 1))
    assert xml.value.category == "RESOURCE_LIMIT"
    payload = '<img src=x onerror="alert(1)">'
    escaped = components.html_escape(payload)
    assert "<img" not in escaped
    assert "&lt;" in escaped


def test_msa_rejects_job_id_urls_and_does_not_use_shell():
    import inspect

    source = inspect.getsource(msa._execute_local)
    assert "shell=False" in source
    assert "shell=True" not in source
    with pytest.raises(msa.MsaError) as exc:
        msa.poll_msa(
            {
                "backend": "ebi_clustalo",
                "job_id": "../etc/passwd",
                "submitted_monotonic": 0.0,
            },
            "user@example.com",
            now=10_000,
        )
    assert exc.value.category == "INVALID_INPUT"


def test_rnafold_subprocess_is_allowlisted_and_does_not_use_shell():
    import inspect

    source = inspect.getsource(rna_folding._fold_executable)
    assert "shell=False" in source
    assert "shell=True" not in source
    assert "--noPS" in source
    assert "--noconv" in source
    assert "--temp=" in source
    if not rna_folding.tool_availability()["available"]:
        with pytest.raises(rna_folding.FoldingError) as missing:
            rna_folding.fold_rna("ACGUACGU")
        assert missing.value.category == "TOOL_NOT_INSTALLED"


def test_casoffinder_subprocess_is_allowlisted_and_does_not_use_shell():
    import inspect

    source = inspect.getsource(crispr_casoffinder._run_queries)
    assert "shell=False" in source
    assert "shell=True" not in source
    detect = inspect.getsource(crispr_casoffinder.detect_cas_offinder)
    assert "CAS_OFFINDER_NAMES" in detect or "cas-offinder" in detect


def test_protein_structure_rejects_user_and_loopback_urls():
    with pytest.raises(protein_structure.StructureError) as loopback:
        protein_structure.fetch_structure_file("https://127.0.0.1/secret.cif")
    assert loopback.value.category == "INVALID_INPUT"
    with pytest.raises(protein_structure.StructureError) as user_url:
        protein_structure.fetch_structure_file("https://evil.example/1CRN.cif")
    assert user_url.value.category == "INVALID_INPUT"


def test_phylogeny_subprocess_and_taxonomy_query_are_constrained():
    import inspect

    fast = inspect.getsource(phylogeny._run_fasttree)
    iq = inspect.getsource(phylogeny._run_iqtree)
    assert "shell=False" in fast and "shell=True" not in fast
    assert "shell=False" in iq and "shell=True" not in iq
    with pytest.raises(ncbi_fetch.NCBIQueryError) as url:
        ncbi_fetch.fetch_taxonomy("https://127.0.0.1/taxonomy", "user@example.com")
    assert url.value.category == "invalid_input"
    with pytest.raises(phylogeny.PhylogenyError) as huge:
        phylogeny.parse_newick("A" * (phylogeny.MAX_NEWICK_CHARS + 1))
    assert huge.value.category == "RESOURCE_LIMIT"
    with pytest.raises(phylogeny.PhylogenyError) as bad:
        phylogeny.parse_newick("((((broken")
    assert bad.value.category == "PARSING_ERROR"
    assert protein_structure.url_is_allowed("http://files.rcsb.org/download/1CRN.cif") is False


def test_gff_seqid_and_cas_offinder_args_are_constrained():
    from modules import genome_annotation, crispr_casoffinder, opencl_runtime

    with pytest.raises(genome_annotation.AnnotationError):
        genome_annotation.parse_gff_text(
            "..\\windows\ttestdb\tgene\t1\t2\t.\t+\t.\tID=x\n",
            declared_assembly="GRCh38.p14",
        )
    env = opencl_runtime.subprocess_environ(
        {"OCL_ICD_FILENAMES": "Intel_OpenCL_ICD64.dll"}
    )
    assert not opencl_runtime.ocl_icd_filenames_is_hanging(str(env.get("OCL_ICD_FILENAMES") or ""))
    policy = crispr_casoffinder.cas_offinder_version_policy()
    assert policy["v3_production_ready"] is False
