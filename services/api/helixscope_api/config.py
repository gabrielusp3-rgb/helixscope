"""API/system configuration. Scientific defaults stay in helixscope_core."""

from __future__ import annotations

import os

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict

_VERCEL_ORIGIN_KEYS: tuple[str, ...] = (
    "VERCEL_URL",
    "VERCEL_BRANCH_URL",
    "VERCEL_PROJECT_PRODUCTION_URL",
)


def _vercel_preview_origins() -> list[str]:
    """HTTPS origins injected by Vercel. Never '*'. Empty when not on Vercel."""
    origins: list[str] = []
    for key in _VERCEL_ORIGIN_KEYS:
        host = (os.environ.get(key) or "").strip().rstrip("/")
        if not host:
            continue
        origin = host if host.startswith(("https://", "http://")) else f"https://{host}"
        if origin not in origins:
            origins.append(origin)
    return origins


class ApiSettings(BaseSettings):
    """Non-secret API process settings. NCBI keys stay in the environment unused here."""

    model_config = SettingsConfigDict(
        env_prefix="HELIXSCOPE_API_",
        extra="ignore",
    )

    environment: str = "development"
    cors_origins: str = "http://localhost:3000"
    max_upload_bytes: int = 1_048_576
    max_sequence_chars: int = 250_000
    job_max_concurrent: int = 2
    job_queue_capacity: int = 8
    job_default_timeout_s: float = 120.0
    log_level: str = "INFO"

    def cors_origin_list(self) -> list[str]:
        """Parse comma-separated origins plus Vercel deployment hosts. Empty is deny-all (no *)."""
        items = [part.strip() for part in self.cors_origins.split(",")]
        origins = [item for item in items if item and item != "*"]
        for origin in _vercel_preview_origins():
            if origin not in origins:
                origins.append(origin)
        return origins
