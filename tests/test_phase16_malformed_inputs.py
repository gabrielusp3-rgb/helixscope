"""Matriz de entradas malformadas: FASTA, VCF, GFF, mmCIF e Newick.

Cada parser precisa falhar de forma explicita e categorizada. Uma entrada
corrompida que atravessa silenciosamente e pior do que um erro: produz um
resultado cientifico com aparencia normal e origem invalida.

Sem rede. Sem subprocesso. Ficheiros temporarios do pytest.
"""

from __future__ import annotations

import pytest

from modules import (
    genome_annotation,
    genome_fasta,
    phylogeny,
    protein_structure,
    variant_core,
)

MALFORMED_FAI = (
    ("chr1\t100\t6\t60", "too few fields"),
    ("chr1\t100\t6\t60\t61\t9", "too many fields"),
    ("chr1\tNOTANUMBER\t6\t60\t61", "non-integer length"),
    ("chr1\t100\t6\t0\t61", "zero linebases"),
    ("chr1\t100\t6\t60\t10", "linewidth below linebases"),
    ("chr1\t-100\t6\t60\t61", "negative length"),
    ("chr1\t100\t6\t60\t61\nchr1\t200\t120\t60\t61", "duplicate contig name"),
    ("\t100\t6\t60\t61", "empty contig name"),
)


@pytest.mark.parametrize("text,label", MALFORMED_FAI)
def test_malformed_fai_is_rejected_with_a_parsing_error(text, label):
    with pytest.raises(genome_fasta.GenomeFastaError) as error:
        genome_fasta.parse_fai_text(text)
    assert error.value.category == "PARSING_ERROR", label


def test_valid_fai_still_parses_and_keeps_contig_geometry():
    parsed = genome_fasta.parse_fai_text("chr1\t100\t6\t60\t61\nchrM\t16569\t120\t70\t71")
    assert parsed["n_contigs"] == 2
    assert parsed["contigs"]["chrM"]["length"] == 16569
    assert parsed["contigs"]["chr1"]["linebases"] == 60


def test_fasta_with_sequence_before_any_header_is_rejected(tmp_path):
    fasta = tmp_path / "headerless.fa"
    fasta.write_text("ACGTACGTAC\n>chr1\nACGT\n", encoding="ascii")
    with pytest.raises(genome_fasta.GenomeFastaError) as error:
        genome_fasta.build_fai_for_fasta(str(fasta), fai_path=str(tmp_path / "out.fai"))
    assert error.value.category in {"PARSING_ERROR", "INVALID_INPUT"}


def test_fasta_with_ragged_line_lengths_is_rejected(tmp_path):
    fasta = tmp_path / "ragged.fa"
    fasta.write_text(">chr1\nACGTACGTAC\nACG\nACGTACGTAC\n", encoding="ascii")
    with pytest.raises(genome_fasta.GenomeFastaError) as error:
        genome_fasta.build_fai_for_fasta(str(fasta), fai_path=str(tmp_path / "ragged.fai"))
    assert error.value.category == "PARSING_ERROR"


def test_missing_fasta_path_is_invalid_input_not_a_crash(tmp_path):
    with pytest.raises(genome_fasta.GenomeFastaError) as error:
        genome_fasta.build_fai_for_fasta(str(tmp_path / "absent.fa"))
    assert error.value.category == "INVALID_INPUT"


MALFORMED_VCF_LIKE = (
    "17 43093557 C",
    "17 43093557",
    "17 C G",
    "17 notaposition C G",
    "17 43093557 C G extra field",
    "17 4.5 C G",
    "17 1e9 C G",
    "17 43093557 XYZ G",
    "17 43093557 C XYZ",
    "17 43093557 . .",
    "17\t43093557\tC",
)


@pytest.mark.parametrize("text", MALFORMED_VCF_LIKE)
def test_malformed_vcf_like_variant_lines_are_rejected(text):
    with pytest.raises(variant_core.VariantInputError):
        variant_core.build_variant(text=text, assembly="GRCh38.p14")


MALFORMED_HGVS = (
    "NC_000017.11:g.C>G",
    "NC_000017.11:g.43093557",
    "NC_000017.11:43093557C>G",
    "g.43093557C>G",
    "NC_000017.11:g.43093557C>",
    "NC_000017.11:g.0C>G",
)


@pytest.mark.parametrize("text", MALFORMED_HGVS)
def test_malformed_genomic_hgvs_is_rejected_locally(text):
    with pytest.raises(variant_core.VariantInputError):
        variant_core.build_variant(text=text, assembly="GRCh38.p14")


MALFORMED_SPDI = (
    "NC_000017.11:43093556:C",
    "NC_000017.11:notaposition:C:G",
    "NC_000017.11:-1:C:G",
    "NC_000017.11:43093556:C:C",
)


@pytest.mark.parametrize("text", MALFORMED_SPDI)
def test_malformed_spdi_is_rejected(text):
    with pytest.raises(variant_core.VariantInputError):
        variant_core.build_variant(text=text, assembly="GRCh38.p14")


MALFORMED_GFF = (
    ("chr1\tsrc\tgene\t1\n", "too few columns"),
    ("chr1\tsrc\tgene\tnotanumber\t500\t.\t+\t.\tID=g1\n", "non-integer start"),
    ("chr1\tsrc\tgene\t500\t100\t.\t+\t.\tID=g1\n", "end before start"),
    ("chr1\tsrc\tgene\t0\t500\t.\t+\t.\tID=g1\n", "zero start in a 1-based format"),
    ("chr1\tsrc\tgene\t1\t500\t.\tsideways\t.\tID=g1\n", "invalid strand"),
    ("chr1\tsrc\tCDS\t1\t500\t.\t+\t7\tID=c1\n", "invalid phase"),
)


@pytest.mark.parametrize("text,label", MALFORMED_GFF)
def test_malformed_gff_records_are_rejected(text, label):
    with pytest.raises(genome_annotation.AnnotationError) as error:
        genome_annotation.parse_gff_text(text, declared_assembly="GRCh38.p14")
    assert error.value.category in {"PARSING_ERROR", "INVALID_INPUT"}, label


def test_gff_seqid_with_path_or_shell_characters_is_rejected():
    for seqid in ("../../etc/passwd", "chr1; rm -rf /", "chr1\x00"):
        record = f"{seqid}\tsrc\tgene\t1\t500\t.\t+\t.\tID=g1\n"
        with pytest.raises(genome_annotation.AnnotationError):
            genome_annotation.parse_gff_text(record, declared_assembly="GRCh38.p14")


MALFORMED_NEWICK = (
    "(A,B",
    "A,B);",
    "((A,B),C));",
    "",
    "   ",
    "(A,B)C;extra",
)


@pytest.mark.parametrize("text", MALFORMED_NEWICK)
def test_malformed_newick_is_rejected(text):
    with pytest.raises((phylogeny.PhylogenyError, ValueError)):
        phylogeny.parse_newick(text)


def test_valid_newick_still_parses_with_named_leaves():
    tree = phylogeny.parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
    labels = {leaf.name for leaf in tree.get_terminals()}
    assert labels == {"A", "B", "C"}


def test_unnamed_newick_leaves_stay_unnamed_instead_of_getting_invented_labels():
    """Newick sem rotulos e sintaticamente valido mas cientificamente inutil."""
    tree = phylogeny.parse_newick("(:1.0,:2.0);")
    names = [leaf.name for leaf in tree.get_terminals()]
    assert len(names) == 2
    assert all(name in (None, "") for name in names)


MALFORMED_MMCIF = (
    "",
    "data_TEST\n",
    "data_TEST\nloop_\n_atom_site.group_PDB\nATOM\n",
    "data_TEST\nloop_\n_atom_site.Cartn_x\nNaN\n",
    "not a structure file at all",
)


@pytest.mark.parametrize("text", MALFORMED_MMCIF)
def test_malformed_mmcif_is_rejected_instead_of_producing_an_empty_structure(text):
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_mmcif(text)


def test_mmcif_with_nonfinite_coordinates_is_rejected():
    text = (
        "data_TEST\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 N N THR A 1 NaN 14.0 5.0 1 A 1\n"
    )
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_mmcif(text)


def test_malformed_pdb_without_atom_records_is_rejected():
    with pytest.raises(protein_structure.StructureError):
        protein_structure.parse_pdb("HEADER    NOTHING HERE\nEND\n")
