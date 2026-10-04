# Scientific methods

Statuses below are the ones HelixScope stores. They are defined in
[PROVENANCE.md](PROVENANCE.md).

## DNA

Local, in `modules/dna_analysis.py`.

| Method | What it is | What it is not |
| --- | --- | --- |
| Base composition | Counts of A, C, G, T after whitespace is removed and letters are uppercased | A quality-trimmed read set |
| GC and AT content | Percent of G+C or A+T among canonical bases. Ambiguous bases are excluded from the denominator. No canonical base yields NaN, not zero | A melting experiment |
| GC and AT skew | (G-C)/(G+C) and (A-T)/(A+T), including windowed series | A called origin of replication |
| Cumulative skew | Lobry sum: +1 for G, -1 for C, and the same idea for A and T. The final value uses every base. The plot may be downsampled | The skew ratio |
| Reverse complement | Watson-Crick complement, reversed | A translation |
| Melting temperature | Wallace-style and SantaLucia nearest-neighbour reports, with the salt and oligo concentrations shown on the result | A measured Tm |
| Molecular weight | Sum of residue masses for DNA or RNA | An intact-polymer mass from mass spectrometry |
| k-mers | Counts for k of 1 to 5 | A genome-wide k-mer index |
| ORFs | ATG to stop, minimum 150 nt when used as a CDS candidate, with redundant overlaps removed | A gene annotation. Official CDS are preferred when a GenBank accession is fetched |
| CpG | CpG island scan and CpG dinucleotide positions, including observed/expected for that scan | A methylation measurement |
| Shannon entropy | Information content of the residue frequencies | A complexity mask such as SEG |
| Restriction sites | Matches against the embedded enzyme table | A commercial cloning validation |

Input longer than 250,000 raw characters, or 200,000 residues after cleaning,
is rejected. The container checks used a 150,000 nt ACGT repeat for linear
statistics. That length is a validation fixture, not a second biological law.
Plots are capped at 2,500 points. The illustrative DNA helix is capped at
400 nt and is not the analytical sequence.

## RNA

Composition uses the DNA module where the alphabet allows U.

Secondary structure uses ViennaRNA 2.7.2 in the container:

- Global MFE fold is limited to 600 nt. Above that, HelixScope raises
  `RESOURCE_LIMIT` and does not return an empty structure as if it were a fold.
- Ensemble metrics (centroid, diversity, ensemble free energy) use the
  partition function for sequences inside that limit. Status is `PREDICTED`.
- Local analysis calls ViennaRNA `pfl_fold` / `pfl_fold_up` with a default
  window of 70. It does not call the global fold. The record sets
  `global_structure` to false. The technical cap for this local path is
  150,000 nt.

RNA secondary-structure prediction is not an experimental 3D structure and
not an RNA fold from a crystal or cryo-EM model.

## Protein

`modules/protein_analysis.py` uses Biopython ProtParam for standard amino acids.

| Method | Implementation | Status |
| --- | --- | --- |
| Molecular weight | Average residue masses, peptide bonds accounted for. Reported in daltons and kilodaltons | `COMPUTED` |
| Theoretical pI | Bjellqvist pKa scale via ProtParam | `COMPUTED`, not an experimental pI |
| Extinction at 280 nm | ProtParam, reduced cysteines and cystine bridges reported separately | `COMPUTED` |
| Net charge | Henderson-Hasselbalch via ProtParam at the pH you set | `COMPUTED` |
| Shannon entropy | Residue-frequency entropy in bits. This is not SEG | `COMPUTED` |
| GRAVY, instability index, aliphatic index, aromaticity | The corresponding ProtParam or published scale, named on the result | `COMPUTED` |

HelixScope does not predict a protein fold from sequence. A 3D view requires
coordinates from a file or a remote entry.

## Alignment

`modules/alignment.py`.

- Global: Needleman-Wunsch with Gotoh affine gaps.
- Local: Smith-Waterman with the same gap treatment.
- DNA and RNA: match +2, mismatch -1, gap open -5, gap extend -0.5.
- Protein: BLOSUM62.

Each sequence is limited to 8,000 residues. Identity, score and gap counts
come from that alignment. A pairwise alignment is not BLAST and not an MSA.
The dot plot is limited to 2,000 residues per sequence.

## MSA

`modules/msa.py`.

The backend that is always offered is EMBL-EBI Clustal Omega (`ebi_clustalo`).
Local Clustal Omega, MAFFT or MUSCLE appear only when those executables are
installed. They are not installed in the tested container.

Caps: 30 sequences, 8,000 residues each, 80,000 residues in total, 300 s
wait. Offline, the EBI job cannot be submitted. A pre-aligned FASTA can still
be used as input to phylogeny; that file is not a new MSA computed offline.

An MSA is not a tree.

## Phylogeny

`modules/phylogeny.py`, from a completed alignment.

| Method | Where it runs | Notes |
| --- | --- | --- |
| Neighbor-joining | Local | Distance models include p-distance, Jukes-Cantor 1969 and Kimura 1980 for nucleotides |
| UPGMA | Local | Assumes a molecular clock. The drawing root is not a confirmed biological root |
| FastTree | FastTree 2.1.11 in the container | Approximate maximum likelihood. Absent binary returns unavailable, not a fake tree |
| IQ-TREE | IQ-TREE 3.1.3 in the container | Refuses p-distance. Jukes-Cantor and Kimura 1980 are accepted for the tested nucleotide fixture |

Bootstrap for NJ and UPGMA is optional and capped at 50 replicates. Zero
replicates means no support values. IQ-TREE UFBoot or SH-aLRT appear only when
that run requested them and the executable produced them. HelixScope does not
invent branch support.

A tree is a statement about the alignment, the method and the model. It is
not a taxonomic authority tree. Taxonomy from NCBI is not this inference.

## BLAST

Remote BLAST uses the NCBI BLAST URL API and returns NCBI hit data.
Local BLAST+ 2.17.0+ searches the tiny nucleotide and protein fixtures built
into the image (`helixscope_tiny_nucl`, `helixscope_tiny_prot`). Those
fixtures are not nt or nr.

The product `blastn` path uses BLAST+ DUST. An ACGT repeat can be `NO_HITS`
even when `tiny_ref` is the same repeat. That is the low-complexity filter,
not a parser failure. A command-line run with `-dust no` on that fixture is a
separate check; it is not the default product behavior.

BLAST similarity is not CRISPR specificity and not an MSA.

## Structure

- mmCIF and PDB text are parsed locally. Coordinates that fail the parser are
  a parsing error, not a repaired model.
- The 3D view draws those coordinates. It does not predict them.
- US-align 20260328 aligns two monomer files. RMSD and the two TM-scores are
  the numbers US-align prints. The 1CRN-versus-1CRN fixture scored RMSD 0.0 Å
  and TM-score 1.0. That checks a self-alignment. It does not certify other
  pairs.
- RCSB Alignment API results are retrieved.
- AlphaFold files are predicted.
- The short DNA/RNA helix drawing is illustrative and is capped at 400 nt.

## CRISPR

`modules/crispr.py` builds guide tables for the Cas systems declared there,
including SpCas9 with PAM NGG on the 3' side of a 20 nt spacer. Other systems
in that table (SaCas9, Cas12a, Cas12b, Cas13, base editors, prime editor) use
the PAM and the metrics that module defines for them.

The on-target figure is a positional heuristic. The code states that it is
not Doench Rule Set 2 and not DeepHF. Genome-wide CFD and MIT specificity are
not computed. A search of a pasted sequence is a search of that sequence.

Cas-OFFinder is not installed in the tested container. Full-genome search
stays blocked. HelixScope does not simulate hits.

## Variant

`modules/variant_core.py` builds an identity from an assembly and a coordinate
or an HGVS-style token. GRCh37 and GRCh38 coordinates are not interchangeable.
The identity is `COMPUTED` as a coordinate object. It is not a phenotype.

VEP, ClinVar, InterPro and UniProt layers are `RETRIEVED` when the call
succeeds and `UNAVAILABLE` or `UNMAPPED` when it does not. None of those
layers is a HelixScope pathogenicity score.
