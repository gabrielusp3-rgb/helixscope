# Continuous integration

These workflows run on GitHub for every matching push and pull request. A
green run is the record of that commit. The workflows do not deploy this
repository and do not publish a container image.

## Workflows

| Workflow | When it runs | What a failure means |
| --- | --- | --- |
| CI | Every push and pull request | The product pytest suite failed, or a README link, `CITATION.cff` key, or workflow YAML check failed |
| CodeQL | Push, pull request, and Monday 06:00 UTC | The Python analysis uploaded a blocking alert under the repository's code scanning policy |
| Security | Push, pull request, and Monday 06:00 UTC | `pip-audit` found a known issue in `requirements.txt`, or a pull request adds a critical dependency advisory |
| Docker | Changes to the Dockerfile, Compose, requirements, or `docker/`, plus Monday 07:00 UTC and manual dispatch | The image did not build, the pinned Python/Streamlit/Biopython versions drifted, or Trivy reported a critical OS finding |

Skips inside pytest are skips. A missing FastTree, IQ-TREE, BLAST+, US-align, or ViennaRNA installation on the GitHub runner does not become a pass. Live engine tests skip and print the reason, including the ViennaRNA Python checks for local pairing, ensemble metrics, and the circular RNA layout. The overview check still runs; it only requires the FastTree or ViennaRNA Python label when that tool is actually detected. The Docker workflow builds the image that contains those engines and checks Python 3.14.3, Streamlit 1.58.0, Biopython 1.87, ViennaRNA 2.7.2, and the FastTree, IQ-TREE, BLAST+, and US-align binaries. It does not run the full scientific suite. That suite stays a local container run.

Browser end-to-end checks stay on a local machine. AppTest cases that already live under `tests/` run as part of pytest. They are not a browser suite.

## What CI installs

The product job installs `requirements.txt` only. It does not install FastAPI and it does not collect `tests/api`.

`pip-audit` is installed in the Security workflow. It is not added to the Docker image.

## Dependency updates

Dependabot opens weekly pull requests for pip, Docker, and GitHub Actions. It does not merge them. Streamlit, Biopython, and base-image updates still have to pass CI and a scientific check before anyone merges them. The frozen Next.js tree is not a Dependabot ecosystem here.

## Actions

External actions are pinned to a full commit SHA from the official repository, with the release tag in a comment:

- `actions/checkout` v7.0.1
- `actions/setup-python` v7.0.0
- `github/codeql-action` v4.38.2
- `actions/dependency-review-action` v5.0.0
- `aquasecurity/trivy-action` v0.36.0

Workflows request `contents: read`. CodeQL also requests `security-events: write` so it can upload results. There is no `pull_request_target` and no repository secret in these files.

## Repository settings

These are settings, not files in this tree:

- secret scanning and push protection
- branch protection: pull request required, this CI required, no force-push
- code scanning visible from CodeQL

`CODEOWNERS` was not added. There is no GitHub team name to put in it.

## Local equivalents

```text
python -m pytest tests --ignore=tests/api -q
python -m pip install PyYAML
python scripts/check_public_docs.py
docker compose config --quiet
```

A full `docker compose build` is the same command the Docker workflow runs. It is slow because it compiles and downloads the scientific engines.
