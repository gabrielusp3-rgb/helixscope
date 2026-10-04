"""Request-id, body limits, and conservative security headers. Sequences are never logged."""

from __future__ import annotations

import logging
import time
import uuid

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

from helixscope_api.errors import error_body

LOGGER = logging.getLogger("helixscope_api.access")


class RequestIdMiddleware(BaseHTTPMiddleware):
    """Assign an opaque UUID4 request id. Does not hash sequence bodies."""

    async def dispatch(self, request: Request, call_next) -> Response:
        incoming = request.headers.get("x-request-id", "").strip()
        request_id = incoming if incoming and len(incoming) <= 64 else str(uuid.uuid4())
        request.state.request_id = request_id
        started = time.perf_counter()
        response = await call_next(request)
        duration_ms = int((time.perf_counter() - started) * 1000)
        LOGGER.info(
            "request_id=%s method=%s path=%s status=%s duration_ms=%s",
            request_id,
            request.method,
            request.url.path,
            response.status_code,
            duration_ms,
        )
        response.headers["X-Request-ID"] = request_id
        return response


class BodyLimitMiddleware(BaseHTTPMiddleware):
    """Reject oversized Content-Length. Does not trust the header as the only check for uploads."""

    def __init__(self, app, max_bytes: int) -> None:
        super().__init__(app)
        self.max_bytes = int(max_bytes)

    async def dispatch(self, request: Request, call_next) -> Response:
        raw = request.headers.get("content-length")
        if raw:
            try:
                length = int(raw)
            except ValueError:
                length = -1
            if length > self.max_bytes:
                rid = str(getattr(request.state, "request_id", "") or "")
                headers = {"X-Request-ID": rid} if rid else {}
                return JSONResponse(
                    status_code=413,
                    content=error_body("RESOURCE_LIMIT", "Request body exceeds the configured size limit.", rid),
                    headers=headers,
                )
        return await call_next(request)


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Conservative headers. No CSP theater; no credentialed wildcard CORS."""

    async def dispatch(self, request: Request, call_next) -> Response:
        response = await call_next(request)
        response.headers.setdefault("X-Content-Type-Options", "nosniff")
        response.headers.setdefault("X-Frame-Options", "DENY")
        response.headers.setdefault("Referrer-Policy", "no-referrer")
        response.headers.setdefault("Cache-Control", "no-store")
        return response
