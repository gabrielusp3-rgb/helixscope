"""Scientific status vocabulary. Strings match ``modules.provenance``.

Do not mix with job execution states (QUEUED, RUNNING, ...).
"""

from __future__ import annotations

from modules.provenance import EVIDENCE_STATUSES, STRUCTURE_KINDS

SCIENTIFIC_STATUSES: tuple[str, ...] = EVIDENCE_STATUSES

__all__ = ("SCIENTIFIC_STATUSES", "STRUCTURE_KINDS", "EVIDENCE_STATUSES")
