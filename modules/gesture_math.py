"""Camera-orbit math for optional Hand Control. No scientific coordinates.

Hand landmarks drive yaw, pitch and radius of a Plotly scene camera.
This module never mutates atom, residue or alignment coordinates.

Implementacao heuristica simplificada inspirada em Casiez et al. 2012 One Euro
Filter; nao reproduz o modelo/algoritmo original publicado.
"""

from __future__ import annotations

import math
from typing import Mapping, Optional, Sequence

WRIST: int = 0
THUMB_TIP: int = 4
INDEX_TIP: int = 8
INDEX_MCP: int = 5
PINKY_MCP: int = 17
MIDDLE_MCP: int = 9

PINCH_ENGAGE: float = 0.38
PINCH_RELEASE: float = 0.50
DEAD_ZONE: float = 0.012
YAW_GAIN: float = 2.8
PITCH_GAIN: float = 2.4
PITCH_MIN: float = -1.15
PITCH_MAX: float = 1.15
RADIUS_MIN: float = 0.55
RADIUS_MAX: float = 7.5
ZOOM_GAIN: float = 1.8
MAX_YAW_STEP: float = 0.18
MAX_PITCH_STEP: float = 0.14
MIN_PALM: float = 1e-4
MIN_CONFIDENCE: float = 0.55
RESET_HOLD_S: float = 2.0
OPEN_PALM_PINCH: float = 0.85
OPEN_PALM_MIN_CONFIDENCE: float = 0.75
OPEN_PALM_STABILITY: float = 0.035
GESTURE_RESET_DEFAULT: bool = False

Landmark = Mapping[str, float]


def _finite(value: object, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(number):
        return default
    return number


def landmark_xyz(points: Sequence[Landmark], index: int) -> tuple[float, float, float]:
    """Return one landmark as (x, y, z). Missing points are (0, 0, 0).

    Args:
        points: MediaPipe-style landmarks with x/y/z.
        index: Landmark index.

    Returns:
        Cartesian triple. Missing index is zeros, not an exception.

    Raises:
        Nenhum.
    """
    if index < 0 or index >= len(points):
        return (0.0, 0.0, 0.0)
    item = points[index] or {}
    return (_finite(item.get("x")), _finite(item.get("y")), _finite(item.get("z")))


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def palm_scale(points: Sequence[Landmark]) -> float:
    """Hand-size normalizer: index MCP to pinky MCP, else wrist to middle MCP.

    Args:
        points: One hand, 21 landmarks when complete.

    Returns:
        Positive scale. Falls back to 1.0 if degenerate.

    Raises:
        Nenhum.
    """
    scale = _dist(landmark_xyz(points, INDEX_MCP), landmark_xyz(points, PINKY_MCP))
    if scale < MIN_PALM:
        scale = _dist(landmark_xyz(points, WRIST), landmark_xyz(points, MIDDLE_MCP))
    return scale if scale >= MIN_PALM else 1.0


def pinch_ratio(points: Sequence[Landmark]) -> float:
    """Thumb-index distance divided by palm scale.

    Args:
        points: One hand.

    Returns:
        Dimensionless ratio. Smaller is more pinched.

    Raises:
        Nenhum.
    """
    return _dist(landmark_xyz(points, THUMB_TIP), landmark_xyz(points, INDEX_TIP)) / palm_scale(
        points
    )


def tracking_ok(confidence: object) -> bool:
    """Whether a hand detection is stable enough to drive the camera.

    Args:
        confidence: Detector score in 0..1.

    Returns:
        True when confidence meets MIN_CONFIDENCE.

    Raises:
        Nenhum.
    """
    return _finite(confidence) >= MIN_CONFIDENCE


def is_open_palm(points: Sequence[Landmark], *, grabbed: bool) -> bool:
    """Open-hand candidate for a held reset. False while pinched.

    Args:
        points: One hand.
        grabbed: Current pinch-grab state.

    Returns:
        True when the hand is open and not grabbed.

    Raises:
        Nenhum.
    """
    if grabbed:
        return False
    return pinch_ratio(points) >= OPEN_PALM_PINCH


def palm_centroid(points: Sequence[Landmark]) -> tuple[float, float, float]:
    """Mean of wrist, index MCP and pinky MCP. Used for hold stability.

    Args:
        points: One hand.

    Returns:
        Cartesian centroid.

    Raises:
        Nenhum.
    """
    wrist = landmark_xyz(points, WRIST)
    index = landmark_xyz(points, INDEX_MCP)
    pinky = landmark_xyz(points, PINKY_MCP)
    return (
        (wrist[0] + index[0] + pinky[0]) / 3.0,
        (wrist[1] + index[1] + pinky[1]) / 3.0,
        (wrist[2] + index[2] + pinky[2]) / 3.0,
    )


def open_palm_reset_ready(
    *,
    enabled: bool,
    grabbed: bool,
    points: Sequence[Landmark],
    confidence: object,
    elapsed_s: float,
    start_centroid: Optional[tuple[float, float, float]],
) -> bool:
    """Whether an open-palm hold should fire camera reset.

    Args:
        enabled: Explicit UI opt-in. Default off.
        grabbed: Pinch-grab active.
        points: Current landmarks.
        confidence: Detector score.
        elapsed_s: Seconds the open-palm candidate has been held.
        start_centroid: Palm centroid when the hold started.

    Returns:
        True only when opt-in, high confidence, open hand, stable, and held.

    Raises:
        Nenhum.
    """
    if not enabled or grabbed:
        return False
    if _finite(confidence) < OPEN_PALM_MIN_CONFIDENCE:
        return False
    if not is_open_palm(points, grabbed=grabbed):
        return False
    if elapsed_s < RESET_HOLD_S:
        return False
    if start_centroid is None:
        return False
    now = palm_centroid(points)
    displacement = _dist(now, start_centroid) / palm_scale(points)
    if displacement > OPEN_PALM_STABILITY:
        return False
    return True


def apply_pinch_hysteresis(ratio: float, grabbed: bool) -> bool:
    """Engage/release pinch with separate thresholds.

    Args:
        ratio: Current pinch_ratio.
        grabbed: Previous grab state.

    Returns:
        New grab state.

    Raises:
        Nenhum.
    """
    if grabbed:
        return ratio <= PINCH_RELEASE
    return ratio <= PINCH_ENGAGE


def two_hand_span(
    left: Sequence[Landmark],
    right: Sequence[Landmark],
) -> float:
    """Normalized distance between wrists, divided by mean palm scale.

    Args:
        left: First hand.
        right: Second hand.

    Returns:
        Positive span. Degenerate hands yield 1.0.

    Raises:
        Nenhum.
    """
    scale = 0.5 * (palm_scale(left) + palm_scale(right))
    if scale < MIN_PALM:
        scale = 1.0
    return _dist(landmark_xyz(left, WRIST), landmark_xyz(right, WRIST)) / scale


def eye_to_orbit(
    eye: Mapping[str, float],
    center: Optional[Mapping[str, float]] = None,
) -> dict[str, float]:
    """Convert Plotly eye (Z-up) to yaw, pitch, radius about center.

    Args:
        eye: Camera eye x/y/z.
        center: Orbit center. Default origin.

    Returns:
        Dict yaw, pitch, radius. Pitch is clamped.

    Raises:
        Nenhum.
    """
    cx = _finite((center or {}).get("x"))
    cy = _finite((center or {}).get("y"))
    cz = _finite((center or {}).get("z"))
    dx = _finite(eye.get("x")) - cx
    dy = _finite(eye.get("y")) - cy
    dz = _finite(eye.get("z")) - cz
    radius = math.sqrt(dx * dx + dy * dy + dz * dz)
    if radius < 1e-8:
        return {"yaw": 0.0, "pitch": 0.4, "radius": 2.4}
    pitch = math.asin(max(-1.0, min(1.0, dz / radius)))
    yaw = math.atan2(dy, dx)
    return {
        "yaw": yaw,
        "pitch": max(PITCH_MIN, min(PITCH_MAX, pitch)),
        "radius": max(RADIUS_MIN, min(RADIUS_MAX, radius)),
    }


def orbit_to_eye(
    yaw: float,
    pitch: float,
    radius: float,
    center: Optional[Mapping[str, float]] = None,
) -> dict[str, float]:
    """Rebuild Plotly eye from spherical orbit. Z-up. Never returns NaN.

    Args:
        yaw: Horizontal angle radians.
        pitch: Elevation radians, clamped.
        radius: Distance, clamped.
        center: Orbit center.

    Returns:
        Dict x, y, z.

    Raises:
        Nenhum.
    """
    safe_pitch = max(PITCH_MIN, min(PITCH_MAX, _finite(pitch)))
    safe_radius = max(RADIUS_MIN, min(RADIUS_MAX, _finite(radius, 2.4)))
    if not math.isfinite(yaw):
        yaw = 0.0
    cos_p = math.cos(safe_pitch)
    cx = _finite((center or {}).get("x"))
    cy = _finite((center or {}).get("y"))
    cz = _finite((center or {}).get("z"))
    return {
        "x": cx + safe_radius * cos_p * math.cos(yaw),
        "y": cy + safe_radius * cos_p * math.sin(yaw),
        "z": cz + safe_radius * math.sin(safe_pitch),
    }


def apply_dead_zone(dx: float, dy: float, zone: float = DEAD_ZONE) -> tuple[float, float]:
    """Zero small normalized hand deltas.

    Args:
        dx: Horizontal landmark delta.
        dy: Vertical landmark delta.
        zone: Half-width of the dead zone.

    Returns:
        Filtered (dx, dy).

    Raises:
        Nenhum.
    """
    fx = 0.0 if abs(dx) < zone else dx
    fy = 0.0 if abs(dy) < zone else dy
    return (fx, fy)


def clamp_orbit_delta(dyaw: float, dpitch: float) -> tuple[float, float]:
    """Limit per-frame angular steps.

    Args:
        dyaw: Proposed yaw change.
        dpitch: Proposed pitch change.

    Returns:
        Clamped pair.

    Raises:
        Nenhum.
    """
    return (
        max(-MAX_YAW_STEP, min(MAX_YAW_STEP, _finite(dyaw))),
        max(-MAX_PITCH_STEP, min(MAX_PITCH_STEP, _finite(dpitch))),
    )


def orbit_from_pinch_move(
    *,
    yaw: float,
    pitch: float,
    radius: float,
    dx: float,
    dy: float,
    grabbed: bool,
) -> dict[str, float]:
    """Apply one-hand grab motion to the orbit. No-op when not grabbed.

    Args:
        yaw: Current yaw.
        pitch: Current pitch.
        radius: Current radius.
        dx: Normalized hand delta x (right positive).
        dy: Normalized hand delta y (up positive in image space; inverted for pitch).
        grabbed: Whether pinch is engaged.

    Returns:
        Updated yaw, pitch, radius.

    Raises:
        Nenhum.
    """
    if not grabbed:
        return {"yaw": yaw, "pitch": pitch, "radius": radius}
    fx, fy = apply_dead_zone(dx, dy)
    dyaw, dpitch = clamp_orbit_delta(fx * YAW_GAIN, -fy * PITCH_GAIN)
    return {
        "yaw": yaw + dyaw,
        "pitch": max(PITCH_MIN, min(PITCH_MAX, pitch + dpitch)),
        "radius": max(RADIUS_MIN, min(RADIUS_MAX, radius)),
    }


def zoom_from_two_hands(
    *,
    radius_start: float,
    span: float,
    span_ref: float,
) -> float:
    """Map two-hand span to radius. Hands apart -> zoom in (smaller r).

    Call with the radius captured when the two-hand gesture engaged, not the
    previous frame's radius, to avoid compounding.

    Args:
        radius_start: Camera radius at zoom engagement.
        span: Current two_hand_span.
        span_ref: Span at zoom engagement.

    Returns:
        Clamped radius.

    Raises:
        Nenhum.
    """
    if span_ref < MIN_PALM or span < MIN_PALM:
        return max(RADIUS_MIN, min(RADIUS_MAX, radius_start))
    ratio = span / span_ref
    factor = 1.0 / max(0.25, min(4.0, ratio ** ZOOM_GAIN))
    return max(RADIUS_MIN, min(RADIUS_MAX, radius_start * factor))


class OneEuro:
    """Low-latency smoother for noisy landmarks.

    Implementacao heuristica simplificada inspirada em Casiez et al. 2012 One Euro
    Filter; nao reproduz o modelo/algoritmo original publicado.
    """

    def __init__(self, min_cutoff: float = 1.2, beta: float = 0.04, dcutoff: float = 1.0) -> None:
        self.min_cutoff = min_cutoff
        self.beta = beta
        self.dcutoff = dcutoff
        self._x: Optional[float] = None
        self._dx: float = 0.0
        self._t: Optional[float] = None

    def reset(self) -> None:
        """Clear filter memory."""
        self._x = None
        self._dx = 0.0
        self._t = None

    def filter(self, value: float, timestamp_s: float) -> float:
        """Smooth one sample.

        Args:
            value: Raw value.
            timestamp_s: Time in seconds.

        Returns:
            Filtered value. Non-finite input returns the last finite state or 0.

        Raises:
            Nenhum.
        """
        sample = _finite(value)
        if self._x is None or self._t is None:
            self._x = sample
            self._t = timestamp_s
            return sample
        dt = max(1e-3, timestamp_s - self._t)
        self._t = timestamp_s
        dx = (sample - self._x) / dt
        edx = self._dx + _alpha(self.dcutoff, dt) * (dx - self._dx)
        self._dx = edx
        cutoff = self.min_cutoff + self.beta * abs(edx)
        self._x = self._x + _alpha(cutoff, dt) * (sample - self._x)
        return self._x


def _alpha(cutoff: float, dt: float) -> float:
    tau = 1.0 / (2.0 * math.pi * max(1e-4, cutoff))
    return 1.0 / (1.0 + tau / dt)


def traces_xyz_fingerprint(traces: Sequence[Mapping[str, object]]) -> tuple[int, str]:
    """Stable fingerprint of scientific XYZ arrays. Camera is not included.

    Args:
        traces: Plotly-like traces with x, y, z sequences.

    Returns:
        (n_finite_points, comma-joined rounded coordinates).

    Raises:
        Nenhum.
    """
    parts: list[str] = []
    count = 0
    for trace in traces:
        xs = list(trace.get("x") or [])
        ys = list(trace.get("y") or [])
        zs = list(trace.get("z") or [])
        for x, y, z in zip(xs, ys, zs):
            if x is None and y is None and z is None:
                continue
            fx, fy, fz = _finite(x), _finite(y), _finite(z)
            parts.append(f"{fx:.6f}:{fy:.6f}:{fz:.6f}")
            count += 1
    return count, ",".join(parts)


def synthetic_open_hand(*, x: float = 0.5, y: float = 0.5) -> list[dict[str, float]]:
    """Deterministic 21-point open hand for tests. Not a camera capture.

    Args:
        x: Wrist x.
        y: Wrist y.

    Returns:
        21 landmarks.

    Raises:
        Nenhum.
    """
    points = [{"x": x, "y": y, "z": 0.0} for _ in range(21)]
    points[THUMB_TIP] = {"x": x - 0.12, "y": y - 0.02, "z": 0.0}
    points[INDEX_TIP] = {"x": x + 0.02, "y": y - 0.16, "z": 0.0}
    points[INDEX_MCP] = {"x": x + 0.01, "y": y - 0.06, "z": 0.0}
    points[PINKY_MCP] = {"x": x + 0.08, "y": y - 0.05, "z": 0.0}
    points[MIDDLE_MCP] = {"x": x + 0.04, "y": y - 0.06, "z": 0.0}
    return points


def synthetic_pinch_hand(*, x: float = 0.5, y: float = 0.5) -> list[dict[str, float]]:
    """Deterministic pinched hand for tests.

    Args:
        x: Wrist x.
        y: Wrist y.

    Returns:
        21 landmarks with thumb near index.

    Raises:
        Nenhum.
    """
    points = synthetic_open_hand(x=x, y=y)
    index = points[INDEX_TIP]
    points[THUMB_TIP] = {"x": index["x"] + 0.004, "y": index["y"] + 0.004, "z": 0.0}
    return points
