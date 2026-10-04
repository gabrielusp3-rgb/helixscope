"""FastAPI TestClient with lifespan."""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from helixscope_api.main import create_app


@pytest.fixture
def client() -> TestClient:
    with TestClient(create_app()) as test_client:
        yield test_client
