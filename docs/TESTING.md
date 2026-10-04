# Testing

## Product suite

```text
python -m pytest tests --ignore=tests/api -q
```

`tests/api` belongs to the retired FastAPI service. The Streamlit product
suite does not install FastAPI to make that directory pass.

A clean clone of the published `main` tree, without the optional engines,
reported 1200 passed, 25 skipped and 0 failed on 2026-10-04. Earlier runs
inside the tested container, which includes ViennaRNA and the local engines,
reported 1216 passed with 8 skipped and 1215 passed with 9 skipped, both with
0 failed. Those container counts were not repeated on the published commit.
An optional remote test can skip when the network call is unavailable. The
counts describe those runs. They are not a coverage percentage.

The suite includes unit tests, property tests, AppTest checks for navigation
and history, and the container port parser. AppTest is an in-process Streamlit
test. It is not a browser.

## Browser

Playwright, using the Edge channel on this workstation, has been used against
the local container at 1440×900, 1366×768 and 1920×1080. Those checks covered
Overview, DNA, Protein, RNA, the back and forward controls, and the history
panel. They are not committed as a CI job.

## Container checks

Separate from pytest, the running image has been checked for:

- health on ports 8501, 9000 and 10000
- rejection of malformed `PORT` values
- a short scientific script: DNA through 150,000 nt, RNA fold and the 600 nt
  limit, protein indices, pairwise alignment, NJ, UPGMA, FastTree, IQ-TREE,
  local BLAST, 1CRN parsing and US-align
- a 100-iteration memory loop that stayed flat at the sampled RSS
- 2, 4 and 8 browser sessions, with different sequences kept in different sessions

## Security checks in the suite

`tests/test_docker_closure.py` and the phase security tests cover path
payloads, shell metacharacters in parsers, rejected private and link-local
URLs, malformed FASTA/PDB/mmCIF, HTML escaping, and redaction of a fake API
key in job logs. They are not a penetration test of a public deployment.

## What is not claimed

There is no measured line-coverage number in this repository. Do not read the
pass count as 100% coverage.
