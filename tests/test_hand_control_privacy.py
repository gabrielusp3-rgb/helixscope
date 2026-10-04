"""Hand Control privacy and trust-boundary checks. No camera is opened."""

from __future__ import annotations

from pathlib import Path

from ui import hand_control

ROOT = Path(__file__).resolve().parents[1]


def test_hand_control_js_requests_video_only_and_pins_mediapipe() -> None:
    js = hand_control.HAND_CONTROL_JS
    assert "audio: false" in js
    assert "audio: true" not in js
    assert "getUserMedia" in js
    assert "@mediapipe/tasks-vision@" in js
    assert "MP_VERSION = \"0.10.18\"" in js
    assert js.count("0.10.18") >= 1
    assert not any(token in js for token in ('tasks-vision@latest', "/float16/latest/"))
    assert "/float16/1/hand_landmarker.task" in js
    assert "Plotly.relayout" in js
    assert "scene.camera.eye" in js
    assert "trace.x" not in js
    assert "microphone" not in js.lower() or "Microphone is not requested" in js


def test_hand_control_does_not_accept_sequences() -> None:
    import inspect

    params = inspect.signature(hand_control.render_hand_control).parameters
    assert "sequence" not in params
    assert "fasta" not in params
    assert "dna" not in params
    source = inspect.getsource(hand_control._mount)
    assert "data={" in source
    assert "sequence" not in inspect.signature(hand_control.render_hand_control).parameters
    assert "{sequence" not in source
    assert "fasta" not in source.lower()


def test_hand_control_python_pins_model_and_version() -> None:
    assert hand_control.MEDIAPIPE_TASKS_VISION_VERSION == "0.10.18"
    assert "/latest/" not in hand_control.HAND_LANDMARKER_MODEL
    assert hand_control.HAND_CONTROL_JS == (
        ROOT / "ui" / "hand_control.js"
    ).read_text(encoding="utf-8")


def test_js_constants_match_python_gesture_math() -> None:
    from modules import gesture_math as gm

    js = hand_control.HAND_CONTROL_JS
    assert f"PINCH_ENGAGE = {gm.PINCH_ENGAGE}" in js or "PINCH_ENGAGE = 0.38" in js
    assert "PINCH_RELEASE = 0.5" in js or "PINCH_RELEASE = 0.50" in js
    assert "DEAD_ZONE = 0.012" in js
    assert "RADIUS_MIN = 0.55" in js
    assert "RADIUS_MAX = 7.5" in js
    assert f"RESET_HOLD_S = {gm.RESET_HOLD_S}" in js
    assert "OPEN_PALM_MIN_CONFIDENCE = 0.75" in js
    assert "plotHostKey" in js
    assert "gestureResetEnabled" in js
