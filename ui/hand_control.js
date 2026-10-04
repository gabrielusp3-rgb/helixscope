/**
 * HelixScope Hand Control. Camera frames stay in this browser tab.
 * Gestures update Plotly scene.camera only. Scientific XYZ is never written.
 *
 * Constants must match modules/gesture_math.py.
 * MediaPipe: @mediapipe/tasks-vision@0.10.18 (pinned; not a floating latest tag).
 * Model: Google hand_landmarker float16 /1 (not a floating latest tag).
 * Inference stays on the UI thread at INFER_HZ: Tasks Vision 0.10.18 WASM is
 * not a reliable dedicated-worker runtime in Streamlit Components v2 here.
 * Plot lookup uses Streamlit's documented widget-key CSS class (st-key-*),
 * first the plot host container then the chart key.
 */
const MP_VERSION = "0.10.18";
const MP_WASM =
  "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@" + MP_VERSION + "/wasm";
const MP_MODEL =
  "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task";

const WRIST = 0;
const THUMB_TIP = 4;
const INDEX_TIP = 8;
const INDEX_MCP = 5;
const PINKY_MCP = 17;
const MIDDLE_MCP = 9;
const PINCH_ENGAGE = 0.38;
const PINCH_RELEASE = 0.5;
const DEAD_ZONE = 0.012;
const YAW_GAIN = 2.8;
const PITCH_GAIN = 2.4;
const PITCH_MIN = -1.15;
const PITCH_MAX = 1.15;
const RADIUS_MIN = 0.55;
const RADIUS_MAX = 7.5;
const ZOOM_GAIN = 1.8;
const MAX_YAW_STEP = 0.18;
const MAX_PITCH_STEP = 0.14;
const MIN_PALM = 1e-4;
const MIN_CONFIDENCE = 0.55;
const RESET_HOLD_S = 2.0;
const OPEN_PALM_PINCH = 0.85;
const OPEN_PALM_MIN_CONFIDENCE = 0.75;
const OPEN_PALM_STABILITY = 0.035;
const INFER_HZ = 12;
const MOUSE_MUTE_MS = 450;
const LOST_MS = 220;

function finite(value, fallback) {
  const n = Number(value);
  return Number.isFinite(n) ? n : fallback;
}

function dist(a, b) {
  const dx = a[0] - b[0];
  const dy = a[1] - b[1];
  const dz = a[2] - b[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function xyz(points, index) {
  const item = points[index] || {};
  return [finite(item.x, 0), finite(item.y, 0), finite(item.z, 0)];
}

function palmScale(points) {
  let scale = dist(xyz(points, INDEX_MCP), xyz(points, PINKY_MCP));
  if (scale < MIN_PALM) {
    scale = dist(xyz(points, WRIST), xyz(points, MIDDLE_MCP));
  }
  return scale >= MIN_PALM ? scale : 1;
}

function pinchRatio(points) {
  return dist(xyz(points, THUMB_TIP), xyz(points, INDEX_TIP)) / palmScale(points);
}

function palmCentroidJs(points) {
  const wrist = xyz(points, WRIST);
  const index = xyz(points, INDEX_MCP);
  const pinky = xyz(points, PINKY_MCP);
  return [
    (wrist[0] + index[0] + pinky[0]) / 3,
    (wrist[1] + index[1] + pinky[1]) / 3,
    (wrist[2] + index[2] + pinky[2]) / 3,
  ];
}

function twoHandSpan(left, right) {
  const scale = 0.5 * (palmScale(left) + palmScale(right));
  const safe = scale < MIN_PALM ? 1 : scale;
  return dist(xyz(left, WRIST), xyz(right, WRIST)) / safe;
}

function eyeToOrbit(eye, center) {
  const cx = finite(center && center.x, 0);
  const cy = finite(center && center.y, 0);
  const cz = finite(center && center.z, 0);
  const dx = finite(eye && eye.x, 1.75) - cx;
  const dy = finite(eye && eye.y, 1.55) - cy;
  const dz = finite(eye && eye.z, 1.35) - cz;
  const radius = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (radius < 1e-8) {
    return { yaw: 0, pitch: 0.4, radius: 2.4 };
  }
  const pitch = Math.asin(Math.max(-1, Math.min(1, dz / radius)));
  return {
    yaw: Math.atan2(dy, dx),
    pitch: Math.max(PITCH_MIN, Math.min(PITCH_MAX, pitch)),
    radius: Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radius)),
  };
}

function orbitToEye(yaw, pitch, radius, center) {
  const safePitch = Math.max(PITCH_MIN, Math.min(PITCH_MAX, finite(pitch, 0.4)));
  const safeRadius = Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, finite(radius, 2.4)));
  const safeYaw = Number.isFinite(yaw) ? yaw : 0;
  const cosP = Math.cos(safePitch);
  const cx = finite(center && center.x, 0);
  const cy = finite(center && center.y, 0);
  const cz = finite(center && center.z, 0);
  return {
    x: cx + safeRadius * cosP * Math.cos(safeYaw),
    y: cy + safeRadius * cosP * Math.sin(safeYaw),
    z: cz + safeRadius * Math.sin(safePitch),
  };
}

function applyDeadZone(dx, dy) {
  return [Math.abs(dx) < DEAD_ZONE ? 0 : dx, Math.abs(dy) < DEAD_ZONE ? 0 : dy];
}

function clampOrbitDelta(dyaw, dpitch) {
  return [
    Math.max(-MAX_YAW_STEP, Math.min(MAX_YAW_STEP, finite(dyaw, 0))),
    Math.max(-MAX_PITCH_STEP, Math.min(MAX_PITCH_STEP, finite(dpitch, 0))),
  ];
}

function zoomFromTwoHands(radiusStart, span, spanRef) {
  if (spanRef < MIN_PALM || span < MIN_PALM) {
    return Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radiusStart));
  }
  const ratio = span / spanRef;
  const factor = 1 / Math.max(0.25, Math.min(4, Math.pow(ratio, ZOOM_GAIN)));
  return Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radiusStart * factor));
}

function OneEuro(minCutoff, beta, dCutoff) {
  this.minCutoff = minCutoff;
  this.beta = beta;
  this.dCutoff = dCutoff;
  this.x = null;
  this.dx = 0;
  this.t = null;
}
OneEuro.prototype.reset = function () {
  this.x = null;
  this.dx = 0;
  this.t = null;
};
OneEuro.prototype.filter = function (value, timestampS) {
  const sample = finite(value, 0);
  if (this.x === null || this.t === null) {
    this.x = sample;
    this.t = timestampS;
    return sample;
  }
  const dt = Math.max(1e-3, timestampS - this.t);
  this.t = timestampS;
  const dx = (sample - this.x) / dt;
  const aD = alpha(this.dCutoff, dt);
  this.dx = this.dx + aD * (dx - this.dx);
  const cutoff = this.minCutoff + this.beta * Math.abs(this.dx);
  this.x = this.x + alpha(cutoff, dt) * (sample - this.x);
  return this.x;
};
function alpha(cutoff, dt) {
  const tau = 1 / (2 * Math.PI * Math.max(1e-4, cutoff));
  return 1 / (1 + tau / dt);
}

function safePlotKey(raw) {
  return String(raw || "").replace(/[^A-Za-z0-9_-]/g, "");
}

function findPlot(plotKey, hostKey) {
  const host = safePlotKey(hostKey);
  if (host) {
    const wrap = document.querySelector(".st-key-" + host);
    if (wrap) {
      const nested = wrap.querySelector(".js-plotly-plot");
      if (nested) {
        return nested;
      }
    }
  }
  const key = safePlotKey(plotKey);
  if (!key) {
    return null;
  }
  return document.querySelector(".st-key-" + key + " .js-plotly-plot");
}

function cameraOf(gd) {
  const cam = (((gd || {}).layout || {}).scene || {}).camera || {};
  return {
    eye: Object.assign({ x: 1.75, y: 1.55, z: 1.35 }, cam.eye || {}),
    center: Object.assign({ x: 0, y: 0, z: 0 }, cam.center || {}),
    up: Object.assign({ x: 0, y: 0, z: 1 }, cam.up || {}),
  };
}

function applyCamera(gd, eye, center, up) {
  if (!gd || !window.Plotly || typeof window.Plotly.relayout !== "function") {
    return;
  }
  const payload = {
    "scene.camera.eye": {
      x: finite(eye.x, 1.75),
      y: finite(eye.y, 1.55),
      z: finite(eye.z, 1.35),
    },
    "scene.camera.center": {
      x: finite(center.x, 0),
      y: finite(center.y, 0),
      z: finite(center.z, 0),
    },
    "scene.camera.up": {
      x: finite(up.x, 0),
      y: finite(up.y, 0),
      z: finite(up.z, 1),
    },
  };
  window.Plotly.relayout(gd, payload);
}

function handScore(result, index) {
  const handed = (result && (result.handedness || result.handednesses)) || [];
  const row = handed[index] || [];
  const first = row[0] || {};
  return finite(first.score, 0);
}

export default function (component) {
  const parent = component.parentElement;
  const data = component.data || {};
  let plotKey = safePlotKey(data.plotWidgetKey);
  let plotHostKey = safePlotKey(data.plotHostKey);
  let structureId = String(data.structureId || "");
  let showPreview = Boolean(data.showPreview);
  const illustrative = Boolean(data.illustrative);
  const gestureResetEnabled = Boolean(data.gestureResetEnabled);

  const hud = document.createElement("div");
  hud.className = "hs-hand-hud";
  hud.setAttribute("data-hs-hand-hud", "1");
  const statusEl = document.createElement("div");
  const hintEl = document.createElement("div");
  const privacyEl = document.createElement("div");
  const video = document.createElement("video");
  video.setAttribute("playsinline", "true");
  video.setAttribute("muted", "true");
  video.muted = true;
  video.autoplay = true;
  video.className = "hs-hand-preview";
  const resetBtn = document.createElement("button");
  resetBtn.type = "button";
  resetBtn.textContent = "Reset View";
  resetBtn.style.cssText =
    "margin-top:6px;font-size:11px;background:transparent;color:inherit;border:1px solid currentColor;border-radius:6px;padding:4px 8px;cursor:pointer;";
  hud.appendChild(statusEl);
  hud.appendChild(hintEl);
  hud.appendChild(privacyEl);
  hud.appendChild(resetBtn);
  hud.appendChild(video);
  if (parent) {
    parent.appendChild(hud);
  }

  let status = "REQUESTING CAMERA";
  let detail = "";
  let stream = null;
  let landmarker = null;
  let raf = 0;
  let lastInfer = 0;
  let lastVideoTime = -1;
  let stopped = false;
  let grabbed = false;
  let lastHand = null;
  let lastLostAt = 0;
  let calibrated = false;
  let calibUntil = 0;
  let openPalmSince = 0;
  let openPalmOrigin = null;
  let zoomRef = null;
  let applying = false;
  let muteUntil = 0;
  let boundPlot = null;
  let restCamera = null;
  let orbit = { yaw: 0.7, pitch: 0.5, radius: 2.6 };
  let center = { x: 0, y: 0, z: 0 };
  let up = { x: 0, y: 0, z: 1 };
  const fx = new OneEuro(1.2, 0.04, 1.0);
  const fy = new OneEuro(1.2, 0.04, 1.0);

  function setStatus(next, extra) {
    status = next;
    detail = extra || "";
    const nHands = extra && extra.nHands != null ? extra.nHands : null;
    const quality = extra && extra.quality != null ? extra.quality : null;
    const lines = [];
    lines.push("<strong>Hand control · " + status + "</strong>");
    if (illustrative) {
      lines.push("ILLUSTRATIVE geometry. Camera only; coordinates are not experimental.");
    }
    if (nHands != null) {
      lines.push(nHands + (nHands === 1 ? " hand tracked" : " hands tracked"));
    }
    if (quality != null) {
      lines.push("Tracking quality " + quality);
    }
    if (typeof extra === "string" && extra) {
      lines.push(extra);
    }
    statusEl.innerHTML = lines.join("<br>");
    hintEl.textContent = gestureResetEnabled
      ? "Pinch + move · Rotate   Two hands · Zoom   Release · Stop   Open palm hold 2s · Reset"
      : "Pinch + move · Rotate   Two hands · Zoom   Release · Stop   Use Reset View for camera reset";
    privacyEl.textContent =
      status === "OFF" || status === "UNAVAILABLE" || status === "PERMISSION DENIED"
        ? ""
        : "Camera processing: local. Frames are not uploaded.";
    video.classList.toggle("is-on", Boolean(showPreview && stream));
  }

  function stopTracks() {
    if (stream) {
      stream.getTracks().forEach(function (track) {
        try {
          track.stop();
        } catch (err) {
          /* ignore */
        }
      });
      stream = null;
    }
    video.srcObject = null;
    video.classList.remove("is-on");
  }

  function captureRest(gd) {
    const cam = cameraOf(gd);
    restCamera = {
      eye: { x: cam.eye.x, y: cam.eye.y, z: cam.eye.z },
      center: { x: cam.center.x, y: cam.center.y, z: cam.center.z },
      up: { x: cam.up.x, y: cam.up.y, z: cam.up.z },
    };
    center = restCamera.center;
    up = restCamera.up;
    orbit = eyeToOrbit(restCamera.eye, center);
  }

  function onMouseRelayout() {
    if (applying) {
      return;
    }
    muteUntil = performance.now() + MOUSE_MUTE_MS;
    const gd = findPlot(plotKey, plotHostKey);
    if (gd) {
      const cam = cameraOf(gd);
      center = cam.center;
      up = cam.up;
      orbit = eyeToOrbit(cam.eye, center);
    }
  }

  function bindPlot() {
    const gd = findPlot(plotKey, plotHostKey);
    if (!gd || gd === boundPlot) {
      return gd;
    }
    if (boundPlot && boundPlot.removeListener) {
      try {
        boundPlot.removeListener("plotly_relayout", onMouseRelayout);
      } catch (err) {
        /* ignore */
      }
    }
    boundPlot = gd;
    if (gd.on) {
      gd.on("plotly_relayout", onMouseRelayout);
    }
    if (!restCamera) {
      captureRest(gd);
    }
    return gd;
  }

  function resetView() {
    const gd = bindPlot();
    if (!gd || !restCamera) {
      return;
    }
    orbit = eyeToOrbit(restCamera.eye, restCamera.center);
    center = restCamera.center;
    up = restCamera.up;
    applying = true;
    applyCamera(gd, restCamera.eye, restCamera.center, restCamera.up);
    applying = false;
  }

  resetBtn.addEventListener("click", function (event) {
    event.preventDefault();
    resetView();
  });

  function freezeGesture() {
    grabbed = false;
    lastHand = null;
    zoomRef = null;
    openPalmSince = 0;
    openPalmOrigin = null;
    fx.reset();
    fy.reset();
  }

  function driveCamera(hands, scores, nowS) {
    const gd = bindPlot();
    if (!gd) {
      return;
    }
    if (performance.now() < muteUntil) {
      return;
    }
    const usable = [];
    const usableScores = [];
    for (let i = 0; i < hands.length; i += 1) {
      if (scores[i] >= MIN_CONFIDENCE && hands[i] && hands[i].length >= 21) {
        usable.push(hands[i]);
        usableScores.push(scores[i]);
      }
    }
    if (!usable.length) {
      if (!lastLostAt) {
        lastLostAt = performance.now();
      }
      if (performance.now() - lastLostAt >= LOST_MS) {
        freezeGesture();
        setStatus("TRACKING LOST", "Mouse remains available.");
      }
      return;
    }
    lastLostAt = 0;
    if (!calibrated) {
      if (!calibUntil) {
        calibUntil = performance.now() + 1200;
        setStatus("CALIBRATING", "Keep one hand comfortably in view.");
        return;
      }
      if (performance.now() < calibUntil) {
        setStatus("CALIBRATING", "Hand detected. Move naturally to calibrate range.");
        return;
      }
      calibrated = true;
    }

    const quality = scores.length
      ? Math.round(100 * Math.max.apply(null, scores.filter(Number.isFinite))) + "%"
      : "n/a";

    if (usable.length >= 2) {
      grabbed = false;
      lastHand = null;
      const span = twoHandSpan(usable[0], usable[1]);
      if (!zoomRef) {
        zoomRef = { span: span, radius: orbit.radius };
      }
      orbit.radius = zoomFromTwoHands(zoomRef.radius, span, zoomRef.span);
      const eye = orbitToEye(orbit.yaw, orbit.pitch, orbit.radius, center);
      applying = true;
      applyCamera(gd, eye, center, up);
      applying = false;
      setStatus("ACTIVE", { nHands: 2, quality: quality });
      openPalmSince = 0;
      openPalmOrigin = null;
      return;
    }

    zoomRef = null;
    const hand = usable[0];
    const ratio = pinchRatio(hand);
    grabbed = grabbed ? ratio <= PINCH_RELEASE : ratio <= PINCH_ENGAGE;
    const wx = xyz(hand, WRIST)[0];
    const wy = xyz(hand, WRIST)[1];
    const sx = fx.filter(wx, nowS);
    const sy = fy.filter(wy, nowS);

    if (grabbed) {
      openPalmSince = 0;
      openPalmOrigin = null;
      if (lastHand) {
        const raw = applyDeadZone(sx - lastHand.x, lastHand.y - sy);
        const step = clampOrbitDelta(raw[0] * YAW_GAIN, -raw[1] * PITCH_GAIN);
        orbit.yaw += step[0];
        orbit.pitch = Math.max(PITCH_MIN, Math.min(PITCH_MAX, orbit.pitch + step[1]));
        const eye = orbitToEye(orbit.yaw, orbit.pitch, orbit.radius, center);
        applying = true;
        applyCamera(gd, eye, center, up);
        applying = false;
      }
      lastHand = { x: sx, y: sy };
      setStatus("ACTIVE", { nHands: 1, quality: quality });
      return;
    }

    lastHand = { x: sx, y: sy };
    if (
      gestureResetEnabled &&
      finite(usableScores[0], 0) >= OPEN_PALM_MIN_CONFIDENCE &&
      ratio >= OPEN_PALM_PINCH
    ) {
      const now = palmCentroidJs(hand);
      if (!openPalmSince) {
        openPalmSince = performance.now();
        openPalmOrigin = now;
      } else {
        const scale = palmScale(hand);
        const dx = now[0] - openPalmOrigin[0];
        const dy = now[1] - openPalmOrigin[1];
        const dz = now[2] - openPalmOrigin[2];
        const displacement = Math.sqrt(dx * dx + dy * dy + dz * dz) / (scale || 1);
        if (displacement > OPEN_PALM_STABILITY) {
          openPalmSince = performance.now();
          openPalmOrigin = now;
        } else if (performance.now() - openPalmSince >= RESET_HOLD_S * 1000) {
          resetView();
          openPalmSince = 0;
          openPalmOrigin = null;
        }
      }
    } else {
      openPalmSince = 0;
      openPalmOrigin = null;
    }
    setStatus("ACTIVE", { nHands: 1, quality: quality });
  }

  async function loop() {
    if (stopped) {
      return;
    }
    raf = window.requestAnimationFrame(loop);
    bindPlot();
    if (!landmarker || !stream || video.readyState < 2) {
      return;
    }
    const now = performance.now();
    if (now - lastInfer < 1000 / INFER_HZ) {
      return;
    }
    if (video.currentTime === lastVideoTime) {
      return;
    }
    lastVideoTime = video.currentTime;
    lastInfer = now;
    let result;
    try {
      result = landmarker.detectForVideo(video, now);
    } catch (err) {
      setStatus("ERROR", "Hand tracker failed. Mouse remains available.");
      return;
    }
    const hands = (result && result.landmarks) || [];
    const scores = hands.map(function (_h, i) {
      return handScore(result, i);
    });
    driveCamera(hands, scores, now / 1000);
  }

  async function start() {
    if (!window.isSecureContext) {
      setStatus("UNAVAILABLE", "Camera requires HTTPS or localhost.");
      return;
    }
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
      setStatus("UNAVAILABLE", "This browser does not expose a camera API.");
      return;
    }
    setStatus("REQUESTING CAMERA", "Video only. Microphone is not requested.");
    try {
      stream = await navigator.mediaDevices.getUserMedia({
        audio: false,
        video: {
          facingMode: "user",
          width: { ideal: 640 },
          height: { ideal: 480 },
          frameRate: { ideal: 24, max: 30 },
        },
      });
    } catch (err) {
      const name = (err && err.name) || "";
      if (name === "NotAllowedError" || name === "PermissionDeniedError") {
        setStatus("PERMISSION DENIED", "Camera permission denied. Mouse remains available.");
      } else {
        setStatus("UNAVAILABLE", "Camera could not be opened. Mouse remains available.");
      }
      return;
    }
    video.srcObject = stream;
    try {
      await video.play();
    } catch (err) {
      /* autoplay may wait; tracks are live */
    }
    setStatus("CALIBRATING", "Loading hand model.");
    try {
      const mod = await import(
        "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@" + MP_VERSION
      );
      const vision = await mod.FilesetResolver.forVisionTasks(MP_WASM);
      landmarker = await mod.HandLandmarker.createFromOptions(vision, {
        baseOptions: {
          modelAssetPath: MP_MODEL,
          delegate: "GPU",
        },
        runningMode: "VIDEO",
        numHands: 2,
        minHandDetectionConfidence: 0.55,
        minHandPresenceConfidence: 0.5,
        minTrackingConfidence: 0.5,
      });
    } catch (err) {
      stopTracks();
      setStatus("UNAVAILABLE", "Hand control unavailable. Mouse remains available.");
      return;
    }
    setStatus("CALIBRATING", "Keep one hand comfortably in view.");
    raf = window.requestAnimationFrame(loop);
  }

  setStatus("REQUESTING CAMERA");
  start();

  return function () {
    stopped = true;
    if (raf) {
      window.cancelAnimationFrame(raf);
      raf = 0;
    }
    if (boundPlot && boundPlot.removeListener) {
      try {
        boundPlot.removeListener("plotly_relayout", onMouseRelayout);
      } catch (err) {
        /* ignore */
      }
    }
    boundPlot = null;
    freezeGesture();
    stopTracks();
    if (landmarker && typeof landmarker.close === "function") {
      try {
        landmarker.close();
      } catch (err) {
        /* ignore */
      }
    }
    landmarker = null;
    if (hud && hud.parentNode) {
      hud.parentNode.removeChild(hud);
    }
  };
}
