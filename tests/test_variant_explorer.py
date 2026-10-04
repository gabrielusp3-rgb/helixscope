"""Testes da orquestracao Variant Explorer. Offline, com rede injetada.

Cobrem isolamento de falhas entre camadas, propriedades de residuo sem claim de
estabilidade, mapeamento residuo->coluna de MSA, conservacao e divergencia entre
fontes.
"""

from __future__ import annotations

import json
from urllib.parse import urlparse

import pytest

from modules import (
    clinvar_evidence,
    ensembl_vep,
    protein_domains,
    variant_core,
    variant_explorer,
)


class _FakeHandle:
    def __init__(self, body, url: str) -> None:
        self._body = body.encode("utf-8") if isinstance(body, str) else body
        self._url = url
        self.status = 200

    def read(self) -> bytes:
        return self._body

    def geturl(self) -> str:
        return self._url

    def __enter__(self) -> "_FakeHandle":
        return self

    def __exit__(self, *_args) -> bool:
        return False


VEP_MISSENSE = [
    {
        "input": "17 43093557 . C G",
        "assembly_name": "GRCh38",
        "seq_region_name": "17",
        "start": 43093557,
        "end": 43093557,
        "strand": 1,
        "allele_string": "C/G",
        "most_severe_consequence": "missense_variant",
        "transcript_consequences": [
            {
                "transcript_id": "ENST00000357654",
                "gene_id": "ENSG00000012048",
                "gene_symbol": "BRCA1",
                "biotype": "protein_coding",
                "canonical": 1,
                "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
                "protein_start": 658,
                "protein_end": 658,
                "protein_id": "ENSP00000350283",
                "amino_acids": "M/I",
                "codons": "atG/atC",
                "swissprot": ["P38398.282"],
            },
            {
                "transcript_id": "ENST00000471181",
                "gene_symbol": "BRCA1",
                "biotype": "protein_coding",
                "consequence_terms": ["synonymous_variant"],
                "protein_start": 697,
                "protein_end": 697,
                "protein_id": "ENSP00000418960",
                "amino_acids": "M",
                "codons": "atG/atC",
                "swissprot": ["P38398.282"],
            },
            {
                "transcript_id": "ENST00000352993",
                "gene_symbol": "BRCA1",
                "biotype": "protein_coding",
                "consequence_terms": ["intron_variant"],
            },
        ],
    }
]

VEP_INTERGENIC = [
    {
        "input": "2 30000000 . A G",
        "assembly_name": "GRCh38",
        "seq_region_name": "2",
        "start": 30000000,
        "end": 30000000,
        "allele_string": "A/G",
        "most_severe_consequence": "intergenic_variant",
        "intergenic_consequences": [
            {"impact": "MODIFIER", "consequence_terms": ["intergenic_variant"]}
        ],
    }
]

INTERPRO_PAGE = {
    "count": 1,
    "next": None,
    "results": [
        {
            "metadata": {
                "accession": "IPR031099",
                "name": "BRCA1-associated",
                "source_database": "interpro",
                "type": "family",
                "member_databases": {},
            },
            "proteins": [
                {
                    "accession": "p38398",
                    "protein_length": 1863,
                    "entry_protein_locations": [
                        {"fragments": [{"start": 500, "end": 900, "dc-status": "CONTINUOUS"}]}
                    ],
                }
            ],
        }
    ],
}

UNIPROT_ENTRY = {
    "primaryAccession": "P38398",
    "uniProtkbId": "BRCA1_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "genes": [{"geneName": {"value": "BRCA1"}}],
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "sequence": {"length": 1863, "crc64": "ABC"},
    "uniProtKBCrossReferences": [{"database": "AlphaFoldDB", "id": "P38398"}],
}

CLINVAR_RECORD = {
    "uid": "54425",
    "accession": "VCV000054425",
    "accession_version": "VCV000054425.81",
    "title": "NM_007294.4(BRCA1):c.1974G>C (p.Met658Ile)",
    "obj_type": "single nucleotide variant",
    "genes": [{"symbol": "BRCA1", "GeneID": "672"}],
    "germline_classification": {
        "description": "Conflicting classifications of pathogenicity",
        "review_status": "criteria provided, conflicting classifications",
        "last_evaluated": "2026/01/31 00:00",
        "trait_set": [{"trait_name": "Hereditary cancer-predisposing syndrome", "trait_xrefs": []}],
    },
    "supporting_submissions": {"scv": ["SCV1", "SCV2"], "rcv": ["RCV1"]},
    "variation_set": [
        {
            "variation_name": "NM_007294.4(BRCA1):c.1974G>C",
            "cdna_change": "c.1974G>C",
            "variation_xrefs": [],
            "variation_loc": [
                {
                    "status": "current",
                    "assembly_name": "GRCh38",
                    "chr": "17",
                    "start": "43093557",
                    "stop": "43093557",
                    "assembly_acc_ver": "GCF_000001405.38",
                    "ref": "",
                    "alt": "",
                }
            ],
        }
    ],
}


def _vep_opener(payload):
    def opener(request, timeout=None):
        return _FakeHandle(json.dumps(payload), "https://rest.ensembl.org/vep/x")

    return opener


def _clinvar_opener(uids, records):
    def opener(request, timeout=None):
        url = str(request.full_url)
        if "esearch" in url:
            body = {"esearchresult": {"count": str(len(uids)), "idlist": uids, "querytranslation": "t"}}
        else:
            result = {"uids": uids}
            for record in records:
                result[record["uid"]] = record
            body = {"result": result}
        return _FakeHandle(json.dumps(body), "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/x")

    return opener


def _domain_opener():
    def opener(request, timeout=None):
        url = str(request.full_url)
        host = urlparse(url).hostname or ""
        if host == "rest.uniprot.org":
            return _FakeHandle(json.dumps(UNIPROT_ENTRY), "https://rest.uniprot.org/uniprotkb/P38398.json")
        return _FakeHandle(json.dumps(INTERPRO_PAGE), "https://www.ebi.ac.uk/interpro/api/x")

    return opener


@pytest.fixture(autouse=True)
def _no_backoff_sleep(monkeypatch):
    """Remove o sleep do backoff: os testes exercitam a categoria, nao a espera."""
    for module in (ensembl_vep, clinvar_evidence, protein_domains):
        monkeypatch.setattr(module.time, "sleep", lambda *_a, **_k: None)


@pytest.fixture()
def variant() -> dict:
    return variant_core.build_variant(
        text="17 43093557 C G", assembly="GRCh38.p14", accession="GCF_000001405.40"
    )


def test_layer_result_requires_a_reason_when_not_available() -> None:
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.layer_result("domains", "ERROR")
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.layer_result("not_a_layer", "AVAILABLE")
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.layer_result("domains", "WEIRD", reason="x")
    ok = variant_explorer.layer_result("domains", "UNMAPPED", reason="no accession")
    assert ok["status"] == "UNMAPPED"


def test_full_chain_reaches_every_layer(variant) -> None:
    result = variant_explorer.explore_variant(
        variant=variant,
        email="helixscope@example.org",
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
        clinvar_urlopen_fn=_clinvar_opener(["54425"], [CLINVAR_RECORD]),
        domain_urlopen_fn=_domain_opener(),
    )
    summary = result["status_summary"]
    assert summary["variant_identity"] == "AVAILABLE"
    assert summary["consequences"] == "AVAILABLE"
    assert summary["protein_mapping"] == "AVAILABLE"
    assert summary["residue_properties"] == "AVAILABLE"
    assert summary["domains"] == "AVAILABLE"
    assert summary["protein_reference"] == "AVAILABLE"
    assert summary["clinical_evidence"] == "AVAILABLE"
    assert summary["evolutionary_context"] == "SKIPPED"


def test_transcript_ambiguity_is_preserved_not_collapsed(variant) -> None:
    result = variant_explorer.explore_variant(
        variant=variant,
        enable_clinvar=False,
        enable_domains=False,
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
    )
    protein = result["layers"]["protein_mapping"]["data"]
    assert protein["count"] == 2
    ids = [row["transcript_id"] for row in protein["transcripts"]]
    assert "ENST00000357654" in ids and "ENST00000471181" in ids
    assert protein["transcripts"][0]["canonical"] is True
    assert "none is hidden" in protein["note"]


def test_clinvar_failure_does_not_destroy_structural_layers(variant) -> None:
    def failing(request, timeout=None):
        raise OSError("clinvar is down")

    result = variant_explorer.explore_variant(
        variant=variant,
        email="helixscope@example.org",
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
        clinvar_urlopen_fn=failing,
        domain_urlopen_fn=_domain_opener(),
    )
    assert result["status_summary"]["clinical_evidence"] == "ERROR"
    assert result["status_summary"]["domains"] == "AVAILABLE"
    assert result["status_summary"]["protein_mapping"] == "AVAILABLE"
    assert result["status_summary"]["consequences"] == "AVAILABLE"
    assert "not evidence about the variant" in result["layers"]["clinical_evidence"]["reason"]


def test_interpro_failure_does_not_destroy_vep(variant) -> None:
    def failing(request, timeout=None):
        raise OSError("interpro is down")

    result = variant_explorer.explore_variant(
        variant=variant,
        enable_clinvar=False,
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
        domain_urlopen_fn=failing,
    )
    assert result["status_summary"]["domains"] == "ERROR"
    assert result["status_summary"]["protein_reference"] == "ERROR"
    assert result["status_summary"]["consequences"] == "AVAILABLE"
    assert result["status_summary"]["residue_properties"] == "AVAILABLE"


def test_vep_failure_leaves_variant_identity_intact(variant) -> None:
    def failing(request, timeout=None):
        raise OSError("ensembl is down")

    result = variant_explorer.explore_variant(
        variant=variant,
        enable_clinvar=False,
        enable_domains=False,
        vep_urlopen_fn=failing,
    )
    assert result["status_summary"]["variant_identity"] == "AVAILABLE"
    assert result["status_summary"]["consequences"] == "ERROR"
    assert result["status_summary"]["protein_mapping"] == "UNMAPPED"
    assert result["layers"]["variant_identity"]["data"]["identity_hash"]


def test_intergenic_variant_is_unmapped_not_an_error() -> None:
    variant = variant_core.build_variant(text="2 30000000 A G", assembly="GRCh38.p14")
    result = variant_explorer.explore_variant(
        variant=variant,
        enable_clinvar=False,
        enable_domains=False,
        vep_urlopen_fn=_vep_opener(VEP_INTERGENIC),
    )
    protein = result["layers"]["protein_mapping"]
    assert protein["status"] == "UNMAPPED"
    assert "not a failure" in protein["reason"]
    assert result["layers"]["residue_properties"]["status"] == "UNMAPPED"


def test_disabled_layers_say_unknown_not_absent(variant) -> None:
    result = variant_explorer.explore_variant(
        variant=variant, enable_vep=False, enable_clinvar=False, enable_domains=False
    )
    assert "unknown, not absent" in result["layers"]["consequences"]["reason"]
    assert "not evidence that the variant is benign" in result["layers"]["clinical_evidence"]["reason"]


def test_clinvar_without_email_is_skipped_not_reported_as_no_record(variant) -> None:
    result = variant_explorer.explore_variant(
        variant=variant,
        email="",
        enable_domains=False,
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
    )
    layer = result["layers"]["clinical_evidence"]
    assert layer["status"] == "SKIPPED"
    assert "not evidence that the variant is benign" in layer["reason"]


def test_residue_properties_never_claim_stability() -> None:
    change = variant_explorer.residue_property_change("M", "I")
    assert change["status"] == "COMPUTED"
    assert change["stability_claim"] is None
    assert change["hydropathy_difference"] == pytest.approx(2.6)
    assert change["charge_class_changed"] is False
    assert change["volume_difference_a3"] == pytest.approx(3.8)
    text = change["disclaimer"].lower()
    assert "destabilizing" in text and "does not state" in text


def test_residue_properties_detect_charge_and_category_change() -> None:
    change = variant_explorer.residue_property_change("E", "K")
    assert change["charge_class_reference"] == "negative"
    assert change["charge_class_variant"] == "positive"
    assert change["charge_class_changed"] is True
    assert change["side_chain_category_changed"] is True
    glycine = variant_explorer.residue_property_change("G", "W")
    assert glycine["glycine_involved"] is True
    assert glycine["volume_difference_a3"] == pytest.approx(167.7)


def test_residue_properties_refuse_non_standard_residues() -> None:
    for pair in (("X", "A"), ("A", "*"), ("", "A"), ("U", "A")):
        with pytest.raises(variant_explorer.VariantExplorerError):
            variant_explorer.residue_property_change(*pair)


def test_cache_key_changes_with_assembly_allele_and_version(variant) -> None:
    base = dict(variant=variant, layer="consequences", tool="Ensembl VEP", tool_version="116")
    key = variant_explorer.variant_cache_key(**base)
    assert key == variant_explorer.variant_cache_key(**base)
    other_release = variant_explorer.variant_cache_key(**{**base, "tool_version": "115"})
    assert key != other_release
    other_layer = variant_explorer.variant_cache_key(**{**base, "layer": "domains"})
    assert key != other_layer
    grch37 = variant_core.build_variant(text="17 43093557 C G", assembly="GRCh37")
    assert key != variant_explorer.variant_cache_key(**{**base, "variant": grch37})
    other_alt = variant_core.build_variant(text="17 43093557 C T", assembly="GRCh38.p14")
    assert key != variant_explorer.variant_cache_key(**{**base, "variant": other_alt})
    with_params = variant_explorer.variant_cache_key(**base, parameters={"domains": 1})
    assert key != with_params


def test_cache_key_works_for_identifier_only_variants() -> None:
    identifier = variant_core.build_variant(text="rs80357382", assembly="GRCh38.p14")
    key = variant_explorer.variant_cache_key(
        variant=identifier, layer="consequences", tool="Ensembl VEP", tool_version="116"
    )
    assert len(key) == 64


def test_uniprot_accessions_are_extracted_and_validated() -> None:
    row = {"swissprot": ["P38398.282"], "uniprot_isoform": ["P38398-1"], "trembl": ["garbage"]}
    assert variant_explorer.uniprot_accessions_from_consequence(row) == ["P38398"]
    assert variant_explorer.uniprot_accessions_from_consequence({}) == []


def test_residue_maps_to_the_right_alignment_column() -> None:
    rows = ["MKT--ARG", "MK-TAR-G"]
    first = variant_explorer.map_residue_to_msa_column(
        aligned_rows=rows, row_index=0, residue_1based=4
    )
    assert first["status"] == "AVAILABLE"
    assert first["column_1based"] == 6
    assert first["residue_in_alignment"] == "A"
    second = variant_explorer.map_residue_to_msa_column(
        aligned_rows=rows, row_index=1, residue_1based=4
    )
    assert second["column_1based"] == 5


def test_residue_past_the_end_of_a_row_is_unmapped() -> None:
    outcome = variant_explorer.map_residue_to_msa_column(
        aligned_rows=["MKT--ARG"], row_index=0, residue_1based=99
    )
    assert outcome["status"] == "UNMAPPED"
    assert "past the end" in outcome["reason"]


def test_map_residue_validates_its_inputs() -> None:
    for kwargs in (
        dict(aligned_rows=[], row_index=0, residue_1based=1),
        dict(aligned_rows=["ABC"], row_index=5, residue_1based=1),
        dict(aligned_rows=["ABC"], row_index=0, residue_1based=0),
        dict(aligned_rows=["ABC"], row_index="x", residue_1based=1),
    ):
        with pytest.raises(variant_explorer.VariantExplorerError):
            variant_explorer.map_residue_to_msa_column(**kwargs)


def test_evolutionary_context_reports_distribution_and_conservation() -> None:
    rows = ["MKTAR", "MKTAR", "MKTAR", "MKTIR"]
    labels = ["human", "mouse", "chicken", "zebrafish"]
    context = variant_explorer.evolutionary_context_at_residue(
        aligned_rows=rows,
        labels=labels,
        row_index=0,
        residue_1based=4,
        reference_amino_acid="A",
        variant_amino_acid="I",
    )
    assert context["status"] == "AVAILABLE"
    assert context["column_1based"] == 4
    assert context["evidence_status"] == "COMPUTED"
    counts = {item["symbol"]: item["count"] for item in context["residue_distribution"]}
    assert counts == {"A": 3, "I": 1}
    assert context["variant_allele_present_in_labels"] == ["zebrafish"]
    assert context["variant_allele_present_elsewhere"] is True
    assert 0.0 <= context["conservation_score"] <= 1.0
    assert "not evidence of pathogenicity" in context["disclaimer"]


def test_fully_conserved_column_is_wording_only_never_pathogenic() -> None:
    rows = ["MKTAR"] * 5
    context = variant_explorer.evolutionary_context_at_residue(
        aligned_rows=rows, labels=[f"s{i}" for i in range(5)], row_index=0, residue_1based=4
    )
    assert context["conservation_score"] == pytest.approx(1.0)
    assert context["highly_conserved_position"] is True
    assert context["conservation_wording"] == "highly conserved position among the aligned sequences"
    assert "pathogenic" not in context["conservation_wording"]


def test_reference_residue_mismatch_raises_a_warning_not_a_silent_read() -> None:
    rows = ["MKTAR", "MKTGR"]
    context = variant_explorer.evolutionary_context_at_residue(
        aligned_rows=rows,
        labels=["a", "b"],
        row_index=0,
        residue_1based=4,
        reference_amino_acid="W",
    )
    assert context["warnings"]
    assert "different isoform" in context["warnings"][0]


def test_evolutionary_context_reports_gap_fraction() -> None:
    rows = ["MKTAR", "MKT-R", "MKT-R"]
    context = variant_explorer.evolutionary_context_at_residue(
        aligned_rows=rows, labels=["a", "b", "c"], row_index=0, residue_1based=4
    )
    assert context["gap_fraction"] == pytest.approx(2 / 3)
    assert "Gaps excluded" in context["gap_treatment"]


def test_evolutionary_context_requires_matching_labels() -> None:
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.evolutionary_context_at_residue(
            aligned_rows=["MKTAR"], labels=["a", "b"], row_index=0, residue_1based=1
        )


def test_unmapped_residue_yields_unmapped_evolutionary_context() -> None:
    context = variant_explorer.evolutionary_context_at_residue(
        aligned_rows=["MKTAR"], labels=["a"], row_index=0, residue_1based=99
    )
    assert context["status"] == "UNMAPPED"


def test_cross_source_agreement_reports_shared_symbol() -> None:
    agreement = variant_explorer.cross_source_gene_agreement(
        vep={"transcript_consequences": [{"gene_symbol": "BRCA1"}]},
        clinvar={"records": [{"genes": [{"symbol": "BRCA1"}]}]},
        uniprot={"gene_names": ["BRCA1"]},
    )
    assert agreement["agree"] is True
    assert agreement["shared_symbols"] == ["BRCA1"]
    assert "no silent source priority" in agreement["policy"].lower()


def test_cross_source_disagreement_is_shown_not_resolved() -> None:
    agreement = variant_explorer.cross_source_gene_agreement(
        vep={"transcript_consequences": [{"gene_symbol": "NBR2"}]},
        clinvar={"records": [{"genes": [{"symbol": "BRCA1"}]}]},
    )
    assert agreement["agree"] is False
    assert set(agreement["symbols"]) == {"NBR2", "BRCA1"}
    assert "does not pick a winner" in agreement["note"]


def test_cross_source_with_one_source_is_not_comparable() -> None:
    agreement = variant_explorer.cross_source_gene_agreement(
        vep={"transcript_consequences": [{"gene_symbol": "BRCA1"}]}
    )
    assert agreement["agree"] is None
    assert agreement["comparable"] is False


def test_field_provenance_records_lineage() -> None:
    field = variant_explorer.field_provenance(
        "protein_position",
        value=658,
        source="Ensembl VEP",
        evidence_status="RETRIEVED",
        method="VEP REST",
        version="116",
    )
    assert field["evidence_status"] == "RETRIEVED"
    assert field["version"] == "116"
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.field_provenance("x", value=1, source="s", evidence_status="GUESSED")


def test_explore_variant_rejects_a_non_helixscope_variant() -> None:
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.explore_variant(variant={"status": "MADE_UP"})
    with pytest.raises(variant_explorer.VariantExplorerError):
        variant_explorer.explore_variant(variant="not a mapping")


def test_network_disclosures_are_attached_to_every_run(variant) -> None:
    result = variant_explorer.explore_variant(
        variant=variant, enable_vep=False, enable_clinvar=False, enable_domains=False
    )
    disclosures = result["network_disclosures"]
    assert disclosures["ensembl_vep"]["remote"] is True
    assert disclosures["ncbi_clinvar"]["remote"] is True
    assert disclosures["interpro_uniprot"]["remote"] is True


def test_clinvar_layer_marks_records_that_do_not_match_the_coordinate(variant) -> None:
    far = json.loads(json.dumps(CLINVAR_RECORD))
    far["variation_set"][0]["variation_loc"][0]["start"] = "43045750"
    far["variation_set"][0]["variation_loc"][0]["stop"] = "43045750"
    result = variant_explorer.explore_variant(
        variant=variant,
        email="helixscope@example.org",
        enable_domains=False,
        vep_urlopen_fn=_vep_opener(VEP_MISSENSE),
        clinvar_urlopen_fn=_clinvar_opener(["54425"], [far]),
    )
    data = result["layers"]["clinical_evidence"]["data"]
    assert data["identity_matched_count"] == 0
    assert data["records"][0]["identity_match"] is False
