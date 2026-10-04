"""Deterministic Hand Control math. No camera. No scientific XYZ mutation."""

from __future__ import annotations

import math

from modules import gesture_math as gm


def test_pinch_hysteresis_does_not_flicker() -> None:
    mid = (gm.PINCH_ENGAGE + gm.PINCH_RELEASE) / 2
    grabbed = gm.apply_pinch_hysteresis(mid, False)
    assert grabbed is False
    grabbed = gm.apply_pinch_hysteresis(gm.PINCH_ENGAGE - 0.01, False)
    assert grabbed is True
    grabbed = gm.apply_pinch_hysteresis(mid, True)
    assert grabbed is True
    grabbed = gm.apply_pinch_hysteresis(gm.PINCH_RELEASE + 0.01, True)
    assert grabbed is False


def test_dead_zone_swallows_jitter() -> None:
    dx, dy = gm.apply_dead_zone(0.004, -0.003)
    assert dx == 0.0
    assert dy == 0.0
    dx, dy = gm.apply_dead_zone(0.05, 0.0)
    assert dx == 0.05


def test_open_and_pinch_synthetic_hands() -> None:
    open_hand = gm.synthetic_open_hand()
    pinch = gm.synthetic_pinch_hand()
    assert gm.pinch_ratio(open_hand) > gm.PINCH_RELEASE
    assert gm.pinch_ratio(pinch) < gm.PINCH_ENGAGE
    assert gm.is_open_palm(open_hand, grabbed=False) is True
    assert gm.is_open_palm(pinch, grabbed=False) is False
    assert gm.tracking_ok(0.9) is True
    assert gm.tracking_ok(0.1) is False


def test_open_palm_reset_requires_opt_in_hold_and_stability() -> None:
    open_hand = gm.synthetic_open_hand()
    origin = gm.palm_centroid(open_hand)
    assert gm.GESTURE_RESET_DEFAULT is False
    assert (
        gm.open_palm_reset_ready(
            enabled=False,
            grabbed=False,
            points=open_hand,
            confidence=0.99,
            elapsed_s=3.0,
            start_centroid=origin,
        )
        is False
    )
    assert (
        gm.open_palm_reset_ready(
            enabled=True,
            grabbed=False,
            points=open_hand,
            confidence=0.99,
            elapsed_s=0.2,
            start_centroid=origin,
        )
        is False
    )
    assert (
        gm.open_palm_reset_ready(
            enabled=True,
            grabbed=False,
            points=open_hand,
            confidence=0.2,
            elapsed_s=3.0,
            start_centroid=origin,
        )
        is False
    )
    moved = gm.synthetic_open_hand(x=0.9, y=0.9)
    assert (
        gm.open_palm_reset_ready(
            enabled=True,
            grabbed=False,
            points=moved,
            confidence=0.99,
            elapsed_s=3.0,
            start_centroid=origin,
        )
        is False
    )
    assert (
        gm.open_palm_reset_ready(
            enabled=True,
            grabbed=False,
            points=open_hand,
            confidence=0.99,
            elapsed_s=3.0,
            start_centroid=origin,
        )
        is True
    )


def test_orbit_left_right_up_down_are_bounded() -> None:
    start = gm.eye_to_orbit({"x": 1.75, "y": 1.55, "z": 1.35})
    left = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=-0.08,
        dy=0.0,
        grabbed=True,
    )
    right = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=0.08,
        dy=0.0,
        grabbed=True,
    )
    up = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=0.0,
        dy=0.08,
        grabbed=True,
    )
    down = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=0.0,
        dy=-0.08,
        grabbed=True,
    )
    assert left["yaw"] < start["yaw"]
    assert right["yaw"] > start["yaw"]
    assert up["pitch"] <= start["pitch"]
    assert down["pitch"] >= start["pitch"]
    for state in (left, right, up, down):
        assert gm.PITCH_MIN <= state["pitch"] <= gm.PITCH_MAX
        assert gm.RADIUS_MIN <= state["radius"] <= gm.RADIUS_MAX
        eye = gm.orbit_to_eye(state["yaw"], state["pitch"], state["radius"])
        assert all(math.isfinite(eye[axis]) for axis in ("x", "y", "z"))


def test_release_and_tracking_loss_freeze_orbit() -> None:
    start = gm.eye_to_orbit({"x": 1.75, "y": 1.55, "z": 1.35})
    moved = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=0.2,
        dy=0.2,
        grabbed=False,
    )
    assert moved == start
    lost = gm.orbit_from_pinch_move(
        yaw=1.2,
        pitch=0.3,
        radius=2.0,
        dx=9.0,
        dy=9.0,
        grabbed=False,
    )
    assert lost["yaw"] == 1.2
    assert lost["pitch"] == 0.3


def test_two_hand_zoom_changes_radius_only() -> None:
    left = gm.synthetic_open_hand(x=0.3, y=0.5)
    right = gm.synthetic_open_hand(x=0.7, y=0.5)
    closer = gm.synthetic_open_hand(x=0.45, y=0.5)
    span0 = gm.two_hand_span(left, right)
    span_in = gm.two_hand_span(left, gm.synthetic_open_hand(x=0.9, y=0.5))
    span_out = gm.two_hand_span(left, closer)
    r_in = gm.zoom_from_two_hands(radius_start=2.4, span=span_in, span_ref=span0)
    r_out = gm.zoom_from_two_hands(radius_start=2.4, span=span_out, span_ref=span0)
    assert r_in < 2.4
    assert r_out > 2.4
    assert gm.RADIUS_MIN <= r_in <= gm.RADIUS_MAX
    tiny = gm.zoom_from_two_hands(radius_start=2.4, span=1e-9, span_ref=1e-9)
    assert tiny == 2.4


def test_radius_and_pitch_clamps_and_nan_safety() -> None:
    eye = gm.orbit_to_eye(float("nan"), 99.0, -4.0)
    assert all(math.isfinite(eye[axis]) for axis in ("x", "y", "z"))
    orbit = gm.eye_to_orbit({"x": float("nan"), "y": 0.0, "z": 0.0})
    assert math.isfinite(orbit["yaw"])
    assert math.isfinite(orbit["pitch"])
    assert math.isfinite(orbit["radius"])
    dyaw, dpitch = gm.clamp_orbit_delta(50.0, -50.0)
    assert abs(dyaw) <= gm.MAX_YAW_STEP
    assert abs(dpitch) <= gm.MAX_PITCH_STEP


def test_one_euro_smooths_without_nan() -> None:
    filt = gm.OneEuro()
    prev = filt.filter(0.5, 0.0)
    for i in range(1, 12):
        prev = filt.filter(0.5 + (0.2 if i % 2 else -0.2), i * 0.03)
        assert math.isfinite(prev)
    filt.reset()
    assert filt.filter(1.0, 1.0) == 1.0


def test_camera_orbit_does_not_mutate_scientific_traces() -> None:
    traces = [
        {"x": [1.0, 2.0, None], "y": [0.25, 0.5, None], "z": [-1.0, 3.25, None]},
        {"x": [9.5], "y": [8.25], "z": [7.125]},
    ]
    before = gm.traces_xyz_fingerprint(traces)
    start = gm.eye_to_orbit({"x": 1.75, "y": 1.55, "z": 1.35})
    moved = gm.orbit_from_pinch_move(
        yaw=start["yaw"],
        pitch=start["pitch"],
        radius=start["radius"],
        dx=0.1,
        dy=-0.05,
        grabbed=True,
    )
    zoomed = gm.zoom_from_two_hands(radius_start=moved["radius"], span=2.0, span_ref=1.0)
    gm.orbit_to_eye(moved["yaw"], moved["pitch"], zoomed)
    after = gm.traces_xyz_fingerprint(traces)
    assert after == before
    assert after[0] == 3
    traces[0]["x"][0] = 1.0
    assert gm.traces_xyz_fingerprint(traces) == before
