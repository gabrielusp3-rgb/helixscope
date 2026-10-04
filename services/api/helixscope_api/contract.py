"""Sanitized API/runtime identity. No filesystem paths or secrets."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from fastapi import FastAPI

from helixscope_core import CORE_PACKAGE_VERSION
from helixscope_api import API_CONTRACT_VERSION, PRODUCT_VERSION


def dump_openapi(app: FastAPI) -> str:
    """Canonical OpenAPI dump used for identity hashes.

    Args:
        app: FastAPI application.

    Returns:
        JSON text matching ``tests/api/test_openapi.py`` export format.
    """
    return json.dumps(app.openapi(), indent=2, sort_keys=True) + "\n"


def openapi_sha256(app: FastAPI) -> str:
    """SHA-256 of the running OpenAPI document.

    Args:
        app: FastAPI application.

    Returns:
        Hex digest of the canonical dump.
    """
    return hashlib.sha256(dump_openapi(app).encode("utf-8")).hexdigest()


def runtime_identity(app: FastAPI) -> dict[str, Any]:
    """Product/core/API versions plus OpenAPI hash. No private paths.

    Args:
        app: FastAPI application.

    Returns:
        Sanitized identity mapping for frontend contract checks.
    """
    digest = openapi_sha256(app)
    return {
        "product_version": PRODUCT_VERSION,
        "core_version": CORE_PACKAGE_VERSION,
        "api_version": API_CONTRACT_VERSION,
        "openapi_sha256": digest,
        "runtime_build_identity": (
            f"{PRODUCT_VERSION}+{CORE_PACKAGE_VERSION}+{API_CONTRACT_VERSION}+{digest[:16]}"
        ),
    }
