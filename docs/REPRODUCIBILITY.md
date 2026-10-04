# Reproducibility

Results are reproducible within the declared software and runtime, not across
every future remote database.

## Fixed by the repository

- Application version string `0.24.3-19` in `modules/provenance.py`
- `requirements.txt` pins Streamlit 1.58.0, Biopython 1.87, Plotly 6.9.0,
  pandas 2.3.3, NumPy 2.5.3, pytest 9.1.1 and Hypothesis 6.168.3
- The Dockerfile pins Python 3.14.3 on a recorded Trixie digest and checks
  SHA-256 or MD5 for the FastTree source, the IQ-TREE archive, the BLAST+
  archive and the US-align archive
- ViennaRNA is installed as `ViennaRNA==2.7.2`
- Small scientific fixtures live in `tests/fixtures/`, including `1CRN.cif`,
  `1BNA.cif` and `1RNA.cif`

## Fixed only at image-build time

`apt-get upgrade` installs the Debian security candidates available that day.
Two builds on different days can differ in OpenSSL or glibc even though the
base digest is pinned. That is intentional. Package versions are not frozen a
second time in a `.deb` lock file.

## Not bit-reproducible

- NCBI, Ensembl, RCSB, PDB-REDO, UniProt, InterPro and EMBL-EBI responses
- Any BLAST database outside the tiny fixtures shipped in the image
- Wall-clock timestamps stored as retrieval time
- A phylogeny that asks IQ-TREE for a stochastic search beyond the fixed
  fixture settings used in the tests

## How to check

```text
python -m pytest tests --ignore=tests/api -q
```

Or the same command inside the image, with `tests` mounted and `PYTHONPATH`
set to `/app`. A clean clone of published `main`, without optional engines,
reported 1200 passed, 25 skipped and 0 failed on 2026-10-04. Earlier container
runs reported 1216 passed with 8 skipped, then 1215 passed with 9 skipped,
both with 0 failed. Those container counts were not repeated on the published
commit. Repeat the command after you change code. Do not treat any of these
counts as a permanent property of the project.

A self-alignment of `1CRN.cif` is a fixture check: RMSD 0.0 Å and TM-score 1.0
on that pair only.
