import path from "node:path";
import { defineConfig, devices } from "@playwright/test";

const webRoot = process.cwd();
const repoRoot = path.resolve(webRoot, "../..");
const python = process.env.HELIXSCOPE_PYTHON || process.env.PYTHON || "python";
const reuse = process.env.HELIXSCOPE_E2E_REUSE !== "0";
const pythonPath = [repoRoot, path.join(repoRoot, "services", "api")].join(path.delimiter);

export default defineConfig({
  testDir: "./tests/e2e",
  testIgnore: "vercel-staging-smoke.spec.ts",
  fullyParallel: false,
  workers: 1,
  timeout: 120_000,
  expect: { timeout: 20_000 },
  use: {
    baseURL: "http://localhost:3000",
    browserName: "chromium",
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
    launchOptions: {
      args: ["--use-angle=swiftshader", "--ignore-gpu-blocklist", "--enable-webgl"],
    },
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"] } }],
  webServer: [
    {
      command: `${python} -m uvicorn helixscope_api.main:app --host 127.0.0.1 --port 8000`,
      url: "http://127.0.0.1:8000/health/ready",
      cwd: repoRoot,
      timeout: 120_000,
      reuseExistingServer: reuse,
      env: {
        ...process.env,
        PYTHONPATH: pythonPath,
      },
    },
    {
      command: "npx next start --port 3000 --hostname localhost",
      url: "http://localhost:3000",
      timeout: 180_000,
      reuseExistingServer: reuse,
    },
  ],
});
