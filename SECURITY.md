# Security policy

HelixScope is research software for exploratory sequence, structure, alignment
and genome analysis. It is not a clinical device and it does not claim to be
free of defects.

## Supported versions

Security reports are considered for the current working tree and for the
tested container runtime:

- application version `0.24.3-19` (`modules/provenance.py`)
- Python 3.14.3
- Streamlit 1.58.0
- Biopython 1.87
- the Linux image described by the current `Dockerfile`

Older checkouts are not maintained as separate security branches.

## Reporting a vulnerability

If this repository is published on GitHub, open a private security advisory
on the repository. Do not file a public issue for an unfixed vulnerability,
and do not include exploit details, access tokens, or private sequences in a
public report.

There is no separate security email address for this project.

## What to include

- the version or commit you ran
- whether you used the local Streamlit process or the container
- the smallest input that shows the problem, with any private sequence removed
- the result you observed

## Scope

In scope:

- unexpected execution of user input
- reading or writing files outside the working area
- leakage of environment secrets into logs or the interface
- cross-session disclosure of another user's sequence or history
- a dependency or image issue that is reachable from the HelixScope process

Out of scope:

- results from NCBI, RCSB, EMBL-EBI, ClinVar, Ensembl or other remote services
- absence of an optional local engine (the application is expected to report
  that the engine is not installed)
- theoretical scores that HelixScope already labels as unavailable

## Disclosure

Please allow time to reproduce the report before publishing details. This
project does not offer a paid response window or a bounty. A fix, when one is
made, will be described in `CHANGELOG.md`.
