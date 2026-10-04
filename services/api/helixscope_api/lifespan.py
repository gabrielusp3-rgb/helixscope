"""Application lifespan. No scientific engines, downloads, or remote validation at startup."""

from __future__ import annotations

from contextlib import asynccontextmanager
from typing import AsyncIterator

from fastapi import FastAPI

from helixscope_api.config import ApiSettings
from helixscope_api.jobs import LocalMemoryJobRunner


def create_lifespan(settings: ApiSettings):
    """FastAPI lifespan factory bound to explicit settings."""

    @asynccontextmanager
    async def lifespan(app: FastAPI) -> AsyncIterator[None]:
        app.state.settings = settings
        app.state.jobs = LocalMemoryJobRunner(
            max_concurrent=settings.job_max_concurrent,
            queue_capacity=settings.job_queue_capacity,
        )
        try:
            yield
        finally:
            runner = getattr(app.state, "jobs", None)
            if runner is not None:
                runner.shutdown()

    return lifespan
