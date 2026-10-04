import { defineConfig, devices } from "@playwright/test";

const baseURL = process.env.HELIXSCOPE_STAGING_URL;
if (!baseURL) {
  throw new Error("HELIXSCOPE_STAGING_URL is required for staging Playwright.");
}

export default defineConfig({
  testDir: "./tests/e2e",
  testMatch: "vercel-staging-smoke.spec.ts",
  fullyParallel: false,
  workers: 1,
  timeout: 180_000,
  expect: { timeout: 45_000 },
  use: {
    baseURL,
    browserName: "chromium",
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
    launchOptions: {
      args: ["--use-angle=swiftshader", "--ignore-gpu-blocklist", "--enable-webgl"],
    },
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"] } }],
});
