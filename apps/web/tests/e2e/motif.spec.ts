import { expect, test } from "@playwright/test";
import { API } from "./helpers";

test("motif IUPAC hits match API coordinates", async ({ page, request }) => {
  const sequence = "ACGTACGTACGT";
  const pattern = "ACGT";
  const api = await request.post(`${API}/api/v1/motif/search`, {
    data: { sequence, pattern, molecule: "DNA", include_explanation: true },
  });
  const body = await api.json();
  expect(body.result.n_hits).toBeGreaterThan(0);
  await page.goto("/discovery/motif");
  await page.locator("#motif-seq").fill(sequence);
  await page.locator("#pattern").fill(pattern);
  await page.getByRole("button", { name: "Search" }).click();
  await expect(page.getByTestId("metric-n")).toContainText(String(body.result.n_hits));
  await expect(page.getByTestId("motif-hits")).toContainText(String(body.result.hits[0].start));
  await expect(page.getByTestId("motif-hits")).toContainText(String(body.result.hits[0].end));
  await expect(page.getByTestId("motif-hits")).toContainText(String(body.result.hits[0].match));
  await expect(page.getByTestId("motif-meaning")).toContainText("not proven biological function");
});

test("motif no-hit keeps zero and does not invent function", async ({ page, request }) => {
  const sequence = "ACGTACGTACGT";
  const pattern = "TTTT";
  const api = await request.post(`${API}/api/v1/motif/search`, {
    data: { sequence, pattern, molecule: "DNA", include_explanation: true },
  });
  const body = await api.json();
  expect(body.result.n_hits).toBe(0);
  await page.goto("/discovery/motif");
  await page.locator("#motif-seq").fill(sequence);
  await page.locator("#pattern").fill(pattern);
  await page.getByRole("button", { name: "Search" }).click();
  await expect(page.getByTestId("metric-n")).toContainText("0");
  await expect(page.getByTestId("motif-hits")).toContainText("0 hits.");
});
