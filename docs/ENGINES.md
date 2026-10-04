# Engines

Versions are the ones measured in the tested container. Licenses are not
restated here when this repository does not carry the upstream license file.
See [THIRD_PARTY.md](THIRD_PARTY.md).

## Installed in the tested container

### FastTree 2.1.11

- Role: approximate maximum-likelihood trees from a nucleotide alignment.
- Source used by the image: official `FastTree.c`.
- Status: local executable. A successful run is `COMPUTED` phylogeny, not a
  taxonomic tree.
- If the binary is absent: unavailable. No substitute tree is invented.

### IQ-TREE 3.1.3

- Role: maximum-likelihood trees.
- Source: the IQ-TREE 3 Linux release recorded in the Dockerfile. That file
  records the release as GPL-2.0-or-later.
- Status: local executable. p-distance is refused. Jukes-Cantor and Kimura
  1980 are the nucleotide models used in the container checks.
- If the binary is absent: unavailable.

### BLAST+ 2.17.0+

- Role: local `blastn` and `blastp` against the tiny fixtures built in the
  image.
- Source: NCBI BLAST+ executables.
- Status: local. `NO_HITS` is a real empty hit list. The fixtures are not nt
  or nr.
- If the binary is absent: local BLAST is unavailable. Remote NCBI BLAST is a
  different path.

### US-align 20260328

- Role: monomer structure alignment (`-mol`, `-mm 0`).
- Source: https://zhanggroup.org/US-align/
- Status: local. RMSD and both TM-scores are copied from US-align output.
- If the binary is absent: the local alignment is unavailable. The RCSB
  Alignment API is a different path.

### ViennaRNA 2.7.2

- Role: RNA secondary structure through the Python module, not the RNAfold
  executable.
- Status: `PREDICTED`. Global fold stops at 600 nt with `RESOURCE_LIMIT`.
- If the module is absent: folding is unavailable.

## Not installed in the tested container

| Engine | Status | What HelixScope does instead |
| --- | --- | --- |
| Cas-OFFinder | `not_installed` | Blocks genome-wide search. Does not simulate hits |
| Local DSSP | `not_installed` | PDB-REDO DSSP remains a remote HTTP call |
| STRIDE | `not_installed` | No local STRIDE assignment |
| MSMS | `not_installed` | No molecular surface from MSMS |
| EDTSurf | `not_installed` | No local surface |

OpenCL and a GPU are not part of the container. Cas-OFFinder is not installed
to force that stack on.

## Detection versus validation

Seeing a binary on `PATH` is `detected`. `live_validated` requires a local run
that succeeded. `remote_validated` is a remote check. The Overview screen only
treats live and remote validation as product-visible engine rows for tools
that have been validated; a missing tool is not promoted.
