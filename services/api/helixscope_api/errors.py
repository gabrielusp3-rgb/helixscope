"""HTTP mapping of core/domain errors. No stack traces in bodies."""

from __future__ import annotations

import logging
from typing import Any

from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as HTTPException
from starlette.requests import Request

from helixscope_core.errors import (
    AlignmentApiError,
    BlastError,
    ClinVarError,
    ComparativeError,
    CoordinateError,
    CrisprError,
    DomainError,
    EvidenceError,
    FoldingError,
    GenomeDownloadError,
    GenomeJobError,
    GenomeStoreError,
    MsaError,
    NCBIQueryError,
    OpenCLError,
    PhylogenyError,
    StructureError,
    TaxonomyError,
    ThermodynamicsError,
    USAlignError,
    VariantExplorerError,
    VariantInputError,
    VepError,
)

LOGGER = logging.getLogger("helixscope_api")


class ApiError(Exception):
    """Classified HTTP error that is not a core scientific exception."""

    def __init__(
        self,
        status: int,
        code: str,
        message: str,
        details: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.status = int(status)
        self.code = str(code)
        self.details = details or {}


DOMAIN_EXCEPTIONS: tuple[type[BaseException], ...] = (
    AlignmentApiError,
    BlastError,
    ClinVarError,
    ComparativeError,
    CoordinateError,
    CrisprError,
    DomainError,
    EvidenceError,
    FoldingError,
    GenomeDownloadError,
    GenomeJobError,
    GenomeStoreError,
    MsaError,
    NCBIQueryError,
    OpenCLError,
    PhylogenyError,
    StructureError,
    TaxonomyError,
    ThermodynamicsError,
    USAlignError,
    VariantExplorerError,
    VariantInputError,
    VepError,
    ValueError,
    TypeError,
)


def request_id_of(request: Request) -> str:
    """Opaque request id from middleware state."""
    return str(getattr(request.state, "request_id", "") or "")


def error_body(code: str, message: str, request_id: str, details: dict[str, Any] | None = None) -> dict[str, Any]:
    """Stable error envelope."""
    return {
        "error": {
            "code": code,
            "message": str(message),
            "details": details or {},
            "request_id": request_id,
        }
    }


def _category_of(exc: BaseException) -> str:
    category = str(getattr(exc, "category", "") or getattr(exc, "failure_kind", "") or "")
    return category.strip()


def classify_exception(exc: BaseException) -> tuple[int, str]:
    """Return (http_status, error_code) without leaking internals."""
    if isinstance(exc, HTTPException):
        return int(exc.status_code), "HTTP_ERROR"
    category = _category_of(exc).upper().replace(" ", "_")
    kind = str(getattr(exc, "failure_kind", "") or "").lower()
    text = str(exc).lower()
    mapping: dict[str, tuple[int, str]] = {
        "INVALID_INPUT": (400, "INVALID_INPUT"),
        "INVALID": (400, "INVALID_INPUT"),
        "CONFIGURATION": (400, "INVALID_INPUT"),
        "NOT_FOUND": (404, "NO_RECORD"),
        "WRONG_DATABASE": (400, "INVALID_INPUT"),
        "SOURCE_UNAVAILABLE": (503, "REMOTE_UNAVAILABLE"),
        "TIMEOUT": (504, "TIMEOUT"),
        "RATE_LIMITED": (429, "RATE_LIMITED"),
        "TOOL_NOT_INSTALLED": (503, "ENGINE_NOT_INSTALLED"),
        "UNAVAILABLE": (503, "ENGINE_NOT_INSTALLED"),
        "TOOL_FAILED": (502, "ENGINE_FAILED"),
        "JOB_FAILED": (502, "ENGINE_FAILED"),
        "RESOURCE_LIMIT": (409, "RESOURCE_LIMIT"),
        "PARSING_ERROR": (422, "PARSE_ERROR"),
        "PARSE_ERROR": (422, "PARSE_ERROR"),
        "INVALID_FORMAT": (400, "INVALID_FORMAT"),
        "NO_HITS": (404, "NO_HITS"),
        "METHOD_NOT_APPLICABLE": (400, "METHOD_NOT_APPLICABLE"),
        "MAPPING_UNCERTAIN": (422, "MAPPING_UNCERTAIN"),
        "MAPPING_UNAVAILABLE": (422, "MAPPING_UNAVAILABLE"),
        "INTERNAL": (500, "INTERNAL"),
        "CANCEL_NOT_SUPPORTED": (409, "CANCEL_NOT_SUPPORTED"),
        "JOB_NOT_FOUND": (404, "JOB_NOT_FOUND"),
        "QUEUE_FULL": (409, "RESOURCE_LIMIT"),
    }
    if kind == "timeout" or "timeout" in category.lower():
        return 504, "TIMEOUT"
    if kind == "rate_limited" or "429" in text:
        return 429, "RATE_LIMITED"
    if category in mapping:
        return mapping[category]
    if isinstance(exc, (ValueError, TypeError, VariantInputError, CoordinateError, CrisprError)):
        return 400, "INVALID_INPUT"
    if isinstance(exc, NCBIQueryError) and category.lower() == "not_found":
        return 404, "NO_RECORD"
    return 500, "INTERNAL"


def http_error(status: int, code: str, message: str, request_id: str) -> HTTPException:
    """HTTPException with the API error envelope."""
    return HTTPException(status_code=status, detail=error_body(code, message, request_id))


async def domain_exception_handler(request: Request, exc: BaseException) -> JSONResponse:
    """Map classified core errors. Logs INTERNAL with request id only."""
    status, code = classify_exception(exc)
    rid = request_id_of(request)
    if status >= 500:
        LOGGER.error("request_id=%s code=%s type=%s", rid, code, type(exc).__name__)
        message = "The scientific service encountered an internal error."
    else:
        message = str(exc)
    return JSONResponse(status_code=status, content=error_body(code, message, rid))


async def api_error_handler(request: Request, exc: ApiError) -> JSONResponse:
    """Stable envelope for API-layer classified failures."""
    rid = request_id_of(request)
    return JSONResponse(
        status_code=exc.status,
        content=error_body(exc.code, str(exc), rid, exc.details),
    )


async def unhandled_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Last-resort 500. No traceback in the body."""
    rid = request_id_of(request)
    LOGGER.error("request_id=%s unhandled type=%s", rid, type(exc).__name__)
    return JSONResponse(
        status_code=500,
        content=error_body("INTERNAL", "The scientific service encountered an internal error.", rid),
    )


async def validation_exception_handler(request: Request, exc: RequestValidationError) -> JSONResponse:
    """422 without echoing biological payloads."""
    rid = request_id_of(request)
    details = [
        {"type": item.get("type"), "loc": item.get("loc"), "msg": item.get("msg")}
        for item in exc.errors()
    ]
    return JSONResponse(
        status_code=422,
        content=error_body("INVALID_INPUT", "Request validation failed.", rid, {"errors": details}),
    )


async def http_exception_handler(request: Request, exc: HTTPException) -> JSONResponse:
    """Keep the API error envelope for HTTPException."""
    rid = request_id_of(request)
    detail = exc.detail
    if isinstance(detail, dict) and "error" in detail:
        return JSONResponse(status_code=exc.status_code, content=detail)
    return JSONResponse(
        status_code=int(exc.status_code),
        content=error_body("HTTP_ERROR", str(detail), rid),
    )
