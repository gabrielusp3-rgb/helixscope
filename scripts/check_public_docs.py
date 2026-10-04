"""Check citation metadata, README links, and workflow YAML.

Stdlib plus PyYAML, which CI installs for this script only. The script does
not rewrite files.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

CITATION_KEYS: tuple[str, ...] = (
    "cff-version: 1.2.0",
    "title: HelixScope",
    "version: 0.24.3-19",
    "license: MIT",
    "family-names: Fiusa",
    "given-names: Gabriel Rodrigues Pinto",
)

README_TARGETS: tuple[str, ...] = (
    "docs/SCIENTIFIC_METHODS.md",
    "docs/SCIENTIFIC_INTEGRITY.md",
    "docs/PROVENANCE.md",
    "docs/ENGINES.md",
    "docs/DATA_SOURCES.md",
    "docs/ARCHITECTURE.md",
    "docs/DOCKER.md",
    "docs/REPRODUCIBILITY.md",
    "docs/TESTING.md",
    "docs/LIMITATIONS.md",
    "docs/CONFIGURATION.md",
    "docs/VALIDATION.md",
    "docs/THIRD_PARTY.md",
    "docs/CI.md",
    "SECURITY.md",
    "CITATION.cff",
    "CHANGELOG.md",
    "LICENSE",
)


def check(root: Path) -> list[str]:
    """Return human-readable problems. An empty list means the tree passed.

    Args:
        root: Repository root.

    Returns:
        Problem strings. Empty when the checked files are consistent.
    """
    problems: list[str] = []
    citation = root / "CITATION.cff"
    if not citation.is_file():
        problems.append("CITATION.cff is missing")
    else:
        text = citation.read_text(encoding="utf-8")
        for key in CITATION_KEYS:
            if key not in text:
                problems.append(f"CITATION.cff is missing {key}")
        if "doi:" in text.lower():
            problems.append("CITATION.cff must not invent a DOI")
    readme = root / "README.md"
    if not readme.is_file():
        problems.append("README.md is missing")
    else:
        body = readme.read_text(encoding="utf-8")
        for relative in README_TARGETS:
            if relative not in body:
                problems.append(f"README.md does not link {relative}")
            elif not (root / relative).is_file():
                problems.append(f"README target does not exist: {relative}")
    workflow_dir = root / ".github" / "workflows"
    if not workflow_dir.is_dir():
        problems.append(".github/workflows is missing")
    else:
        for path in sorted(workflow_dir.glob("*.yml")):
            try:
                loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
            except yaml.YAMLError as exc:
                problems.append(f"{path.name} is not valid YAML: {exc}")
                continue
            if not isinstance(loaded, dict) or "jobs" not in loaded:
                problems.append(f"{path.name} has no jobs")
            text = path.read_text(encoding="utf-8")
            if "pull_request_target" in text or "write-all" in text:
                problems.append(f"{path.name} uses a forbidden trigger or permission")
    dependabot = root / ".github" / "dependabot.yml"
    if not dependabot.is_file():
        problems.append(".github/dependabot.yml is missing")
    else:
        try:
            loaded = yaml.safe_load(dependabot.read_text(encoding="utf-8"))
        except yaml.YAMLError as exc:
            problems.append(f"dependabot.yml is not valid YAML: {exc}")
        else:
            if not isinstance(loaded, dict) or loaded.get("version") != 2:
                problems.append("dependabot.yml must be version 2")
    return problems


def main() -> int:
    """Print problems and return 0 only when none were found."""
    problems = check(ROOT)
    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1
    print("public docs and workflow YAML passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
