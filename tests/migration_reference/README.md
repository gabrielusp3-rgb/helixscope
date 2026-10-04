# Migration golden reference

Deterministic scientific fixtures captured during Prompt 0 against HelixScope `0.24.3-19`.

These files are **not** live retrievals. Remote snapshots (Compare/RCSB) are labeled as such.

Regenerate with:

```text
python tests/migration_reference/_freeze.py
```

Do not treat regeneration as license to change formulas.

## Files

| File | Content |
|---|---|
| `dna.json` | Short, mixed, high-GC, low-GC, ambiguous, all-N, 50 kb |
| `rna.json` | Transcription, CAI heuristic, ENC UNAVAILABLE on short CDS, ViennaRNA Python fold |
| `protein.json` | Physicochemical report on BLAST-tiny protein sequence |
| `alignment.json` | Needleman-Wunsch and Smith-Waterman |
| `motif.json` | Exact and IUPAC hits, 0-based half-open |
| `phylogeny.json` | 4-taxon PREALIGNED MSA: NJ, UPGMA, IQ-TREE 3.1.3, FastTree 2.2.0 |
| `crispr.json` | SpCas9 guide on pasted target; model availability; TEST REFERENCE identity |
| `variant.json` | GRCh38 vs GRCh37 identity for `17 43093557 C>G` (no VEP) |
| `structures.json` | SHA-256 of vendored 1BNA/1RNA/1CRN mmCIF |
| `compare.json` | 8HSK vs 8HSF REMOTE_VALIDATED snapshot + fixture parser |
| `evidence.json` | Conflict / no confidence score |
| `blast.json` | BLAST+ 2.17.0+ tiny declared nucleotide DB |
| `parity_classes.json` | EXACT / FLOAT_TOLERANCE / STRUCTURAL_ARRAY / SEMANTIC |
| `function_inventory.json` | Public symbols and coupling classes |

## Parity rules

- Status strings, hashes, integer coordinates, sequences: **EXACT**.
- GC of `NNNN` is **NaN**, encoded as `{"__nonfinite__": "NaN"}`, never `0`.
- Phylogeny: compare leaf sets, method, tool, status. Do **not** require identical Newick across NJ/IQ-TREE/FastTree.
- Compare RMSD: float tolerance 0.05 Å; pair count exact 20.
- Remote fixtures stay fixtures. Live tests stay live.
