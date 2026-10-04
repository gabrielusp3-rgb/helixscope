import { expect, test } from "@playwright/test";

test("command palette navigates modules", async ({ page }) => {
  await page.goto("/overview");
  await page.getByTestId("command-open").click();
  await page.getByTestId("command-input").fill("DNA");
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/analysis\/dna/);
});

test("engines and references render API statuses", async ({ page }) => {
  await page.goto("/overview");
  await expect(page.getByTestId("health-live")).toContainText("alive");
  await expect(page.getByTestId("health-ready")).toContainText("ready");
  await page.goto("/engines");
  await expect(page.locator(".hs-page h1")).toHaveText("Scientific Engines");
  await page.goto("/data/references");
  await expect(page.getByTestId("references-table")).toBeVisible();
});

test("API down shows backend unavailable instead of a blank app", async ({ page }) => {
  await page.route("http://127.0.0.1:8000/**", (route) => route.abort());
  await page.goto("/overview");
  await expect(page.getByTestId("backend-unavailable")).toBeVisible();
  await expect(page.locator(".hs-page h1")).toHaveText("HelixScope");
});
