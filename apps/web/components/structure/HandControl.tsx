"use client";

import { useEffect, useRef, useState } from "react";
import {
  applyPinchHysteresis,
  GESTURE_RESET_DEFAULT,
  openPalmResetReady,
  orbitFromPinchMove,
  palmCentroid,
  pinchRatio,
  trackingOk,
  twoHandSpan,
  zoomFromTwoHands,
  type Landmark,
  type Orbit,
} from "@/lib/gesture/math";
import { ScientificWarning } from "@/components/science/SciencePrimitives";

const INFER_HZ = 12;

export function HandControl({
  orbit,
  onOrbit,
  onReset,
}: {
  orbit: Orbit;
  onOrbit: (orbit: Orbit) => void;
  onReset: () => void;
}) {
  const [enabled, setEnabled] = useState(false);
  const [openPalmReset, setOpenPalmReset] = useState(GESTURE_RESET_DEFAULT);
  const [status, setStatus] = useState("OFF");
  const videoRef = useRef<HTMLVideoElement>(null);
  const streamRef = useRef<MediaStream | null>(null);
  const grabRef = useRef(false);
  const lastPos = useRef<{ x: number; y: number } | null>(null);
  const zoomRef = useRef<{ span: number; radius: number } | null>(null);
  const palmHold = useRef<{ t: number; centroid: [number, number, number] } | null>(null);
  const orbitRef = useRef(orbit);

  useEffect(() => {
    orbitRef.current = orbit;
  }, [orbit]);

  useEffect(() => {
    if (!enabled) {
      streamRef.current?.getTracks().forEach((track) => track.stop());
      streamRef.current = null;
      return;
    }
    let cancelled = false;
    let timer: number | undefined;
    let landmarker: { close: () => void; detectForVideo: (video: HTMLVideoElement, now: number) => { landmarks?: unknown; handedness?: Array<Array<{ score?: number }>> } } | undefined;
    void (async () => {
      try {
        const stream = await navigator.mediaDevices.getUserMedia({
          video: true,
          audio: false,
        });
        if (cancelled) {
          stream.getTracks().forEach((track) => track.stop());
          return;
        }
        streamRef.current = stream;
        if (videoRef.current) {
          videoRef.current.srcObject = stream;
          await videoRef.current.play();
        }
        const vision = await import("@mediapipe/tasks-vision");
        const fileset = await vision.FilesetResolver.forVisionTasks("/mediapipe/wasm");
        landmarker = await vision.HandLandmarker.createFromOptions(fileset, {
          baseOptions: {
            modelAssetPath: "/mediapipe/models/hand_landmarker.task",
          },
          numHands: 2,
          runningMode: "VIDEO",
        });
        setStatus("RUNNING (main thread, not live-camera validated)");
        const tick = () => {
          const video = videoRef.current;
          if (!video || cancelled) return;
          const now = performance.now();
          const result = landmarker?.detectForVideo(video, now);
          if (!result) return;
          const hands = (result.landmarks ?? []) as Landmark[][];
          const conf = result.handedness?.[0]?.[0]?.score;
          if (hands.length === 2 && trackingOk(conf)) {
            const span = twoHandSpan(hands[0], hands[1]);
            if (!zoomRef.current) zoomRef.current = { span, radius: orbitRef.current.radius };
            onOrbit({
              ...orbitRef.current,
              radius: zoomFromTwoHands(zoomRef.current.radius, span, zoomRef.current.span),
            });
            grabRef.current = false;
          } else if (hands.length === 1 && trackingOk(conf)) {
            zoomRef.current = null;
            const points = hands[0];
            const ratio = pinchRatio(points);
            grabRef.current = applyPinchHysteresis(ratio, grabRef.current);
            const wrist = points[0] || { x: 0.5, y: 0.5 };
            const x = Number(wrist.x);
            const y = Number(wrist.y);
            if (grabRef.current && lastPos.current) {
              onOrbit(
                orbitFromPinchMove({
                  ...orbitRef.current,
                  dx: x - lastPos.current.x,
                  dy: -(y - lastPos.current.y),
                  grabbed: true,
                }),
              );
            }
            lastPos.current = { x, y };
            const ready = openPalmResetReady({
              enabled: openPalmReset,
              grabbed: grabRef.current,
              points,
              confidence: conf,
              elapsedS: palmHold.current ? (now - palmHold.current.t) / 1000 : 0,
              startCentroid: palmHold.current?.centroid ?? null,
            });
            if (ready) onReset();
            if (!grabRef.current) {
              palmHold.current ??= { t: now, centroid: palmCentroid(points) };
            } else {
              palmHold.current = null;
            }
          } else {
            grabRef.current = false;
            lastPos.current = null;
            zoomRef.current = null;
          }
          timer = window.setTimeout(tick, 1000 / INFER_HZ);
        };
        tick();
      } catch (error) {
        setStatus(`UNAVAILABLE: ${error instanceof Error ? error.message : "camera or model failed"}`);
        setEnabled(false);
      }
    })();
    return () => {
      cancelled = true;
      if (timer) window.clearTimeout(timer);
      try {
        landmarker?.close();
      } catch {
        /* landmarker may already be closed */
      }
      streamRef.current?.getTracks().forEach((track) => track.stop());
      streamRef.current = null;
    };
  }, [enabled, onOrbit, onReset, openPalmReset]);

  return (
    <div className="hs-panel glass-medium floating-control">
      <h3>Hand Control</h3>
      <p className="hs-lede">
        Camera is off by default. Video only. No microphone. Frames stay in the browser. Status remains
        NOT LIVE CAMERA VALIDATED.
      </p>
      <ScientificWarning>NOT LIVE CAMERA VALIDATED. Synthetic gesture tests are the current evidence.</ScientificWarning>
      <label>
        <input type="checkbox" checked={enabled} onChange={(e) => setEnabled(e.target.checked)} />
        Enable camera Hand Control
      </label>
      <label style={{ display: "block", marginTop: 8 }}>
        <input type="checkbox" checked={openPalmReset} onChange={(e) => setOpenPalmReset(e.target.checked)} />
        Enable open-palm reset (off by default)
      </label>
      <div className="hs-actions">
        <button type="button" className="btn-secondary" onClick={onReset}>
          Reset camera
        </button>
      </div>
      <p data-testid="hand-status">{enabled ? status : "OFF"}</p>
      <video ref={videoRef} muted playsInline style={{ width: 160, display: enabled ? "block" : "none" }} />
    </div>
  );
}
