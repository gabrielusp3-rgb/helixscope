HelixScope test fixtures.

BLAST XML files follow the NCBI BlastOutput schema used by FORMAT_TYPE=XML
on the BLAST Common URL API. They are labeled fixtures for offline tests.
They are not live NCBI responses and must not be treated as such.

MSA FASTA files are synthetic alignments for parser and viewer tests.

RNA fixtures: rna_hsa_let7a.fa is the published hsa-let-7a-5p sequence
(MIMAT0000062). rnafold_output_hairpin.txt tests the RNAfold text parser
and is not a live ViennaRNA result for an arbitrary sequence.

Protein structure:

- 1CRN.cif and rcsb_entry_1CRN.json are public RCSB PDB copies for crambin
  (experimental X-ray, 1.5 A). They are used to test parsing and mapping,
  not as a substitute for a live fetch in production.
- alphafold_prediction_P01308.json is public AlphaFold DB metadata for
  INS_HUMAN. The model is predicted, never experimental.
- uniprot_P01308_structure_xrefs.json is a UniProtKB JSON excerpt with PDB
  and AlphaFoldDB cross-references. NMR entries in that file have no
  resolution; tests must keep resolution as N/A, not 0.
- parser_nmr_altloc.cif and parser_legacy.pdb are parser tests. They are
  not depositions. They exist to cover NMR models, altloc, insertion codes,
  missing residues, and legacy PDB ATOM lines.
- parser_auth_offset.pdb is a parser test whose residue numbers are 10 and
  11 (ALA, LYS), so auth_seq_id is not the polymer index.
- parser_best_effort.cif is a parser test: polymer AGKA with ATOM residues
  numbered 10 and 20, so mapping is BEST_EFFORT, not exact.
- parser_two_chains.cif is a parser test: two polypeptide chains (A=AK, B=GL).
  It is not a deposition. Chains are not concatenated.
- rcsb_search_1CRN.json is a Search API response shape fixture.
- 1BNA.cif and rcsb_entry_1BNA.json are public RCSB PDB copies of the
  Dickerson-Drew B-DNA dodecamer (X-ray, 1.9 A). rcsb_entry_1BNA.json is the
  parsed Data API metadata from GET /rest/v1/core/entry/1BNA (2026-08-27),
  not an invented record.
- 1RNA.cif and rcsb_entry_1RNA.json are public RCSB PDB copies of the
  [U(UA)6A]2 RNA duplex (X-ray, 2.25 A). Same provenance pattern as 1BNA.
- rcsb_polymer_entity_4UN3_{1-4}.json are field subsets of GET
  /rest/v1/core/polymer_entity/4UN3/{n} retrieved 2026-08-27. They are not
  the full API dump and contain no Cartesian coordinates.
- dssp_legacy_min.txt tests the DSSP column parser. It is not a live mkdssp
  run on 1CRN.
- helixscope_test_reference.fa is a synthetic multi-contig FASTA for genome
  index and Cas-OFFinder tests. It is not GRCh38, T2T-CHM13 or GRCm39.
