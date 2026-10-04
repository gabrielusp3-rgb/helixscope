/** Camera-orbit math for Hand Control. Never mutates scientific coordinates. */

export const WRIST = 0;
export const THUMB_TIP = 4;
export const INDEX_TIP = 8;
export const INDEX_MCP = 5;
export const PINKY_MCP = 17;
export const MIDDLE_MCP = 9;

export const PINCH_ENGAGE = 0.38;
export const PINCH_RELEASE = 0.5;
export const DEAD_ZONE = 0.012;
export const YAW_GAIN = 2.8;
export const PITCH_GAIN = 2.4;
export const PITCH_MIN = -1.15;
export const PITCH_MAX = 1.15;
export const RADIUS_MIN = 0.55;
export const RADIUS_MAX = 7.5;
export const ZOOM_GAIN = 1.8;
export const MAX_YAW_STEP = 0.18;
export const MAX_PITCH_STEP = 0.14;
export const MIN_PALM = 1e-4;
export const MIN_CONFIDENCE = 0.55;
export const RESET_HOLD_S = 2.0;
export const OPEN_PALM_PINCH = 0.85;
export const OPEN_PALM_MIN_CONFIDENCE = 0.75;
export const OPEN_PALM_STABILITY = 0.035;
export const GESTURE_RESET_DEFAULT = false;

export type Landmark = { x?: number; y?: number; z?: number };
export type Xyz = { x: number; y: number; z: number };
export type Orbit = { yaw: number; pitch: number; radius: number };

function finite(value: unknown, fallback = 0): number {
  const n = typeof value === "number" ? value : Number(value);
  return Number.isFinite(n) ? n : fallback;
}

export function landmarkXyz(points: Landmark[], index: number): [number, number, number] {
  if (index < 0 || index >= points.length) return [0, 0, 0];
  const item = points[index] || {};
  return [finite(item.x), finite(item.y), finite(item.z)];
}

function dist(a: [number, number, number], b: [number, number, number]): number {
  return Math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2);
}

export function palmScale(points: Landmark[]): number {
  let scale = dist(landmarkXyz(points, INDEX_MCP), landmarkXyz(points, PINKY_MCP));
  if (scale < MIN_PALM) {
    scale = dist(landmarkXyz(points, WRIST), landmarkXyz(points, MIDDLE_MCP));
  }
  return scale >= MIN_PALM ? scale : 1;
}

export function pinchRatio(points: Landmark[]): number {
  return dist(landmarkXyz(points, THUMB_TIP), landmarkXyz(points, INDEX_TIP)) / palmScale(points);
}

export function trackingOk(confidence: unknown): boolean {
  return finite(confidence) >= MIN_CONFIDENCE;
}

export function isOpenPalm(points: Landmark[], grabbed: boolean): boolean {
  if (grabbed) return false;
  return pinchRatio(points) >= OPEN_PALM_PINCH;
}

export function palmCentroid(points: Landmark[]): [number, number, number] {
  const wrist = landmarkXyz(points, WRIST);
  const index = landmarkXyz(points, INDEX_MCP);
  const pinky = landmarkXyz(points, PINKY_MCP);
  return [(wrist[0] + index[0] + pinky[0]) / 3, (wrist[1] + index[1] + pinky[1]) / 3, (wrist[2] + index[2] + pinky[2]) / 3];
}

export function applyPinchHysteresis(ratio: number, grabbed: boolean): boolean {
  return grabbed ? ratio <= PINCH_RELEASE : ratio <= PINCH_ENGAGE;
}

export function twoHandSpan(left: Landmark[], right: Landmark[]): number {
  let scale = 0.5 * (palmScale(left) + palmScale(right));
  if (scale < MIN_PALM) scale = 1;
  return dist(landmarkXyz(left, WRIST), landmarkXyz(right, WRIST)) / scale;
}

export function eyeToOrbit(eye: Xyz, center?: Partial<Xyz>): Orbit {
  const cx = finite(center?.x);
  const cy = finite(center?.y);
  const cz = finite(center?.z);
  const dx = finite(eye.x) - cx;
  const dy = finite(eye.y) - cy;
  const dz = finite(eye.z) - cz;
  const radius = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (radius < 1e-8) return { yaw: 0, pitch: 0.4, radius: 2.4 };
  const pitch = Math.asin(Math.max(-1, Math.min(1, dz / radius)));
  return {
    yaw: Math.atan2(dy, dx),
    pitch: Math.max(PITCH_MIN, Math.min(PITCH_MAX, pitch)),
    radius: Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radius)),
  };
}

export function orbitToEye(yaw: number, pitch: number, radius: number, center?: Partial<Xyz>): Xyz {
  const safePitch = Math.max(PITCH_MIN, Math.min(PITCH_MAX, finite(pitch)));
  const safeRadius = Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, finite(radius, 2.4)));
  const safeYaw = Number.isFinite(yaw) ? yaw : 0;
  const cosP = Math.cos(safePitch);
  const cx = finite(center?.x);
  const cy = finite(center?.y);
  const cz = finite(center?.z);
  return {
    x: cx + safeRadius * cosP * Math.cos(safeYaw),
    y: cy + safeRadius * cosP * Math.sin(safeYaw),
    z: cz + safeRadius * Math.sin(safePitch),
  };
}

export function applyDeadZone(dx: number, dy: number, zone = DEAD_ZONE): [number, number] {
  return [Math.abs(dx) < zone ? 0 : dx, Math.abs(dy) < zone ? 0 : dy];
}

export function clampOrbitDelta(dyaw: number, dpitch: number): [number, number] {
  return [
    Math.max(-MAX_YAW_STEP, Math.min(MAX_YAW_STEP, finite(dyaw))),
    Math.max(-MAX_PITCH_STEP, Math.min(MAX_PITCH_STEP, finite(dpitch))),
  ];
}

export function orbitFromPinchMove(args: {
  yaw: number;
  pitch: number;
  radius: number;
  dx: number;
  dy: number;
  grabbed: boolean;
}): Orbit {
  if (!args.grabbed) return { yaw: args.yaw, pitch: args.pitch, radius: args.radius };
  const [fx, fy] = applyDeadZone(args.dx, args.dy);
  const [dyaw, dpitch] = clampOrbitDelta(fx * YAW_GAIN, -fy * PITCH_GAIN);
  return {
    yaw: args.yaw + dyaw,
    pitch: Math.max(PITCH_MIN, Math.min(PITCH_MAX, args.pitch + dpitch)),
    radius: Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, args.radius)),
  };
}

export function zoomFromTwoHands(radiusStart: number, span: number, spanRef: number): number {
  if (spanRef < MIN_PALM || span < MIN_PALM) {
    return Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radiusStart));
  }
  const ratio = span / spanRef;
  const factor = 1 / Math.max(0.25, Math.min(4, ratio ** ZOOM_GAIN));
  return Math.max(RADIUS_MIN, Math.min(RADIUS_MAX, radiusStart * factor));
}

export function tracesXyzFingerprint(traces: Array<{ x?: unknown[]; y?: unknown[]; z?: unknown[] }>): {
  count: number;
  fingerprint: string;
} {
  const parts: string[] = [];
  let count = 0;
  for (const trace of traces) {
    const xs = trace.x ?? [];
    const ys = trace.y ?? [];
    const zs = trace.z ?? [];
    const n = Math.min(xs.length, ys.length, zs.length);
    for (let i = 0; i < n; i += 1) {
      if (xs[i] == null && ys[i] == null && zs[i] == null) continue;
      parts.push(`${finite(xs[i]).toFixed(6)}:${finite(ys[i]).toFixed(6)}:${finite(zs[i]).toFixed(6)}`);
      count += 1;
    }
  }
  return { count, fingerprint: parts.join(",") };
}

export function syntheticOpenHand(x = 0.5, y = 0.5): Landmark[] {
  const points = Array.from({ length: 21 }, () => ({ x, y, z: 0 }));
  points[THUMB_TIP] = { x: x - 0.12, y: y - 0.02, z: 0 };
  points[INDEX_TIP] = { x: x + 0.02, y: y - 0.16, z: 0 };
  points[INDEX_MCP] = { x: x + 0.01, y: y - 0.06, z: 0 };
  points[PINKY_MCP] = { x: x + 0.08, y: y - 0.05, z: 0 };
  points[MIDDLE_MCP] = { x: x + 0.04, y: y - 0.06, z: 0 };
  return points;
}

export function syntheticPinchHand(x = 0.5, y = 0.5): Landmark[] {
  const points = syntheticOpenHand(x, y);
  const index = points[INDEX_TIP];
  points[THUMB_TIP] = { x: (index.x ?? 0) + 0.004, y: (index.y ?? 0) + 0.004, z: 0 };
  return points;
}

export function openPalmResetReady(args: {
  enabled: boolean;
  grabbed: boolean;
  points: Landmark[];
  confidence: unknown;
  elapsedS: number;
  startCentroid: [number, number, number] | null;
}): boolean {
  if (!args.enabled || args.grabbed) return false;
  if (finite(args.confidence) < OPEN_PALM_MIN_CONFIDENCE) return false;
  if (!isOpenPalm(args.points, args.grabbed)) return false;
  if (args.elapsedS < RESET_HOLD_S) return false;
  if (!args.startCentroid) return false;
  const now = palmCentroid(args.points);
  const displacement = dist(now, args.startCentroid) / palmScale(args.points);
  return displacement <= OPEN_PALM_STABILITY;
}
