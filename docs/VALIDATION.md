# Validation

These are checks of the tested container and fixtures. They are not general
performance or accuracy claims.

## Sequences

- DNA GC on ACGT repeats of 10,000, 50,000, 100,000, 125,000 and 150,000 nt
  was 50.0. The 150,000 nt Lobry cumulative totals were 0 and 0. A short
  reverse complement of the ACGT block stayed ACGT, which is expected for
  that repeat.
- A short melting-temperature call on the first 80 nt of that repeat returned
  71.97 °C under the function's default salt.
- RNA `GGGAAACCC` folded as `PREDICTED`, MFE -1.2 kcal/mol, structure
  `(((...)))`. Ensemble metrics returned `PREDICTED`. A 601 nt string raised
  `RESOURCE_LIMIT`. Local `pfl_fold` on an 800 nt RNA returned `PREDICTED`
  with `global_structure` false.
- Protein `ACDEFGHIKLMNPQRSTVWY`: 2396.0 Da, theoretical pI 6.78, reduced
  extinction 6990, net charge -0.12 at the report pH, Shannon entropy about
  4.32 bits.

## Alignment, MSA and trees

- Needleman-Wunsch on `ACGTACGTACGT` versus `ACGTACGGACGT`: `COMPUTED`,
  identity 91.67%.
- The container offers `ebi_clustalo` and no local MSA binary. Phylogeny
  checks used a four-sequence pre-aligned fixture.
- Neighbor-joining, UPGMA, FastTree and IQ-TREE (Jukes-Cantor) each returned
  `COMPUTED`, four leaves and a Newick string. IQ-TREE with p-distance is
  refused.

## BLAST and structure

- Product `blastn` of the ACGT fixture query: `NO_HITS` (DUST).
- Product `blastp` of the tiny protein fixture query: `READY`, one hit.
- `1CRN.cif` parsed as mmCIF with 327 atoms. US-align of that file against
  itself: `COMPUTED`, RMSD 0.0 Å, TM-score 1.0 on both normalizations.

## Application and container

- Pytest on a clean clone of published `main`, without optional engines:
  1200 passed, 25 skipped, 0 failed on 2026-10-04. Earlier container runs
  reported 1216 passed with 8 skipped, then 1215 passed with 9 skipped, both
  with 0 failed. Those container counts were not repeated on the published
  commit.
- Browser checks at 1440×900, 1366×768 and 1920×1080 covered Overview, DNA,
  Protein, RNA, back, forward and history.
- Health succeeded on `PORT` 8501, 9000 and 10000. Malformed `PORT` values
  exited before Streamlit started.
- A 100-iteration loop of small analyses and history clears stayed at the
  same sampled RSS (147,056 KB on that run).
- Two to eight browser sessions opened, and three concurrent sessions kept
  different sequences in their own fields.
- Cas-OFFinder remained `not_installed`.

Docker Scout on that image reported no critical OS findings and three high
findings without a Debian fix (`gcc-14` runtime libraries and `zlib`). See
[SECURITY.md](../SECURITY.md) and [CHANGELOG.md](../CHANGELOG.md).
