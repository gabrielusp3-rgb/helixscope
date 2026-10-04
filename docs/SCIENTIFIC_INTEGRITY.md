# Scientific integrity

HelixScope is written so that a missing result stays missing.

- Absence is not stored as zero. A GC content with no canonical bases is NaN.
- `None` and empty remote payloads are not filled with a local guess.
- `PREDICTED` is not `EXPERIMENTAL`. An RNA fold and an AlphaFold file are
  predictions.
- `HEURISTIC` is not a published trained model. The CRISPR positional score
  is not Doench Rule Set 2, DeepHF, CFD or MIT.
- An MSA is not a phylogenetic tree. Phylogeny starts only from a completed
  alignment and a method you select.
- A taxonomic name from NCBI is not a phylogeny.
- A BLAST hit is not CRISPR specificity and not an alignment of a whole gene
  family.
- RNA secondary structure is not an experimental RNA 3D structure.
- An illustrative helix is not biological coordinates.
- `UNAVAILABLE`, `not_installed` and `RESOURCE_LIMIT` stay in that state.
  HelixScope does not invent the missing engine output.
- A remote timeout or HTTP error is an error. It is not a new biological record.
- `NO_HITS` means the search reported no hits under its own filters. For local
  `blastn`, DUST can hide a low-complexity literal match.
- A self-alignment of PDB 1CRN (RMSD 0.0 Å, TM-score 1.0) checks that US-align
  and the parser agree on one fixture. It is not a claim about other pairs.
- Genome-wide language is reserved for a finished search of a declared public
  assembly. The tested container cannot do that, because Cas-OFFinder is not
  installed.
