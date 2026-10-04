"""Optional webcam Hand Control for HelixScope Plotly 3D viewports.

Trusted Components v2 script. Does not receive sequences, sequence records,
PDB titles, ClinVar text, or other untrusted scientific strings. Camera frames
are not sent to Python. Gestures change Plotly camera only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import streamlit as st

_JS_PATH = Path(__file__).with_name("hand_control.js")
HAND_CONTROL_JS: str = _JS_PATH.read_text(encoding="utf-8")
HAND_CONTROL_HTML: str = '<div class="hs-hand-root" data-hs-hand="1"></div>'
MEDIAPIPE_TASKS_VISION_VERSION: str = "0.10.18"
HAND_LANDMARKER_MODEL: str = (
    "https://storage.googleapis.com/mediapipe-models/hand_landmarker/"
    "hand_landmarker/float16/1/hand_landmarker.task"
)


def _mount(
    *,
    key: str,
    plot_widget_key: str,
    structure_id: str,
    show_preview: bool,
    illustrative: bool,
    gesture_reset_enabled: bool = False,
) -> Any:
    """Mount the Hand Control script. plot_widget_key is structured data, not HTML.

    Args:
        key: Streamlit widget key for this component instance.
        plot_widget_key: Matching Plotly chart key used to find the canvas.
        structure_id: Opaque structure identity for camera reset, not a sequence.
        show_preview: Local debug preview. Default off at the caller.
        illustrative: When true, HUD keeps the ILLUSTRATIVE disclaimer visible.
        gesture_reset_enabled: Open-palm camera reset. Off unless the user opts in.

    Returns:
        Component result or None when Components v2 cannot mount.

    Raises:
        Nenhum.
    """
    try:
        renderer = st.components.v2.component(
            "hs_hand_control",
            html=HAND_CONTROL_HTML,
            js=HAND_CONTROL_JS,
            isolate_styles=False,
        )
        return renderer(
            key=key,
            data={
                "plotWidgetKey": str(plot_widget_key),
                "plotHostKey": f"{plot_widget_key}_host",
                "structureId": str(structure_id),
                "showPreview": bool(show_preview),
                "illustrative": bool(illustrative),
                "gestureResetEnabled": bool(gesture_reset_enabled),
            },
        )
    except Exception:
        return None


def render_hand_control(
    *,
    plot_key: str,
    structure_id: str,
    illustrative: bool = False,
    control_key: str = "",
) -> None:
    """Render Hand Control UI. Camera stays off until Enable Hand Control.

    Args:
        plot_key: Streamlit key of the Plotly chart this controller may orbit.
        structure_id: Opaque id of the current scientific object (hash/id).
        illustrative: Pass True for ILLUSTRATIVE DNA/RNA geometry.
        control_key: Stable widget-key stem so Reset View does not drop the toggle.

    Returns:
        None.

    Raises:
        Nenhum.
    """
    def _stem(raw: str) -> str:
        return "".join(ch for ch in str(raw) if ch.isalnum() or ch in "_-")

    safe_plot = _stem(plot_key)
    safe_ctrl = _stem(control_key) or safe_plot
    if not safe_plot:
        st.caption("Hand control unavailable")
        return
    enabled = st.toggle(
        "Enable Hand Control",
        value=False,
        key=f"{safe_ctrl}_hand_enable",
        help=(
            "Optional. Requests the camera only after you enable it. "
            "Video only. Hand tracking runs in this browser. "
            "Gestures move the 3D camera, not atomic coordinates."
        ),
    )
    if not enabled:
        st.caption("Hand Control · OFF. Camera is not started.")
        return
    preview = st.checkbox(
        "Tracking preview (local, not saved)",
        value=False,
        key=f"{safe_ctrl}_hand_preview",
        help="Optional local skeleton/video overlay. Default off. Nothing is recorded.",
    )
    gesture_reset = st.checkbox(
        "Allow open-palm camera reset",
        value=False,
        key=f"{safe_ctrl}_hand_palm_reset",
        help=(
            "Off by default. When on, an open palm must stay still for 2 seconds "
            "with high tracking confidence before the camera resets. "
            "Reset View remains available."
        ),
    )
    result = _mount(
        key=f"{safe_ctrl}_hand_js",
        plot_widget_key=safe_plot,
        structure_id=str(structure_id or safe_plot),
        show_preview=bool(preview),
        illustrative=bool(illustrative),
        gesture_reset_enabled=bool(gesture_reset),
    )
    if result is None:
        st.caption("Hand control unavailable. Mouse interaction remains available.")
