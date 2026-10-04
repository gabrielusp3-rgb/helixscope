/**
 * Copy pinned @mediapipe/tasks-vision 0.10.18 WASM into public/.
 * Optionally fetch the pinned Hand Landmarker float16/1 model (Apache-2.0 runtime + Google model terms).
 * Never uses a floating "latest" CDN tag.
 */
import { copyFileSync, mkdirSync, readdirSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { writeFileSync } from "node:fs";

const here = dirname(fileURLToPath(import.meta.url));
const wasmSrc = resolve(here, "../node_modules/@mediapipe/tasks-vision/wasm");
const wasmDest = resolve(here, "../public/mediapipe/wasm");
const modelDestDir = resolve(here, "../public/mediapipe/models");
const licenseDest = resolve(here, "../public/mediapipe/LICENSE.txt");

mkdirSync(wasmDest, { recursive: true });
mkdirSync(modelDestDir, { recursive: true });

for (const name of readdirSync(wasmSrc)) {
  copyFileSync(join(wasmSrc, name), join(wasmDest, name));
}

writeFileSync(
  licenseDest,
  [
    "MediaPipe Tasks Vision WASM files are copied from @mediapipe/tasks-vision@0.10.18",
    "(Apache License 2.0). See node_modules/@mediapipe/tasks-vision.",
    "",
    "Hand Landmarker model (optional, fetched only by this script):",
    "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task",
    "Pinned path /1/ — not /latest/. Camera frames never leave the browser.",
    "",
  ].join("\n"),
  "utf8",
);

const modelUrl =
  "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task";
const modelPath = join(modelDestDir, "hand_landmarker.task");
const fetchModel = process.env.HELIXSCOPE_FETCH_MEDIAPIPE_MODEL !== "0";

if (fetchModel) {
  const response = await fetch(modelUrl);
  if (!response.ok) {
    console.warn(`Hand Landmarker model fetch skipped (${response.status}). WASM copied.`);
    process.exit(0);
  }
  const buf = Buffer.from(await response.arrayBuffer());
  writeFileSync(modelPath, buf);
  console.log(`Wrote ${modelPath} (${buf.length} bytes)`);
}

console.log(`Copied MediaPipe WASM to ${wasmDest}`);
