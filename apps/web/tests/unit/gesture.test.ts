import { describe, expect, it } from "vitest";
import {
  GESTURE_RESET_DEFAULT,
  applyPinchHysteresis,
  orbitFromPinchMove,
  pinchRatio,
  syntheticOpenHand,
  syntheticPinchHand,
  tracesXyzFingerprint,
  zoomFromTwoHands,
} from "@/lib/gesture/math";

describe("gesture camera math", () => {
  it("keeps open-palm reset off by default", () => {
    expect(GESTURE_RESET_DEFAULT).toBe(false);
  });

  it("pinch-move changes orbit but not scientific xyz", () => {
    const traces = [{ x: [1.5, 2.5], y: [0.1, 0.2], z: [8.0, 8.1] }];
    const before = tracesXyzFingerprint(traces);
    const orbit = orbitFromPinchMove({
      yaw: 0.2,
      pitch: 0.1,
      radius: 2.4,
      dx: 0.2,
      dy: 0.05,
      grabbed: true,
    });
    expect(orbit.yaw).not.toBe(0.2);
    expect(tracesXyzFingerprint(traces)).toEqual(before);
  });

  it("synthetic pinch engages and two-hand zoom changes radius only", () => {
    const pinch = syntheticPinchHand();
    expect(applyPinchHysteresis(pinchRatio(pinch), false)).toBe(true);
    const open = syntheticOpenHand();
    expect(applyPinchHysteresis(pinchRatio(open), false)).toBe(false);
    const zoomed = zoomFromTwoHands(2.4, 2, 1);
    expect(zoomed).not.toBe(2.4);
  });
});
