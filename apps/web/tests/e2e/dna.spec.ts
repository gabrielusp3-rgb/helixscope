import { expect, test } from "@playwright/test";
import { API, ORF_150, fillControlled } from "./helpers";

test("DNA analyze matches API values", async ({ page, request }) => {
  const sequence = "ATGC".repeat(30);
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence, include_explanation: true, include_orfs: true, include_profiles: true },
  });
  expect(api.ok()).toBeTruthy();
  const body = await api.json();
  expect(body.result.windowed_profiles.length).toBeGreaterThan(0);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-length")).toContainText(String(body.result.length));
  await expect(page.getByTestId("metric-gc")).toContainText(String(body.result.gc_content));
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
  await expect(page.getByTestId("dna-reverse-complement")).toContainText(String(body.result.reverse_complement));
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
});

test("DNA FASTA paste analyzes residues not the header", async ({ page, request }) => {
  const fasta = ">dna_fixture description with numbers 9606\nATGCATGCATGC\n";
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: fasta, include_explanation: true, include_profiles: true },
  });
  const body = await api.json();
  expect(body.result.status).toBe("COMPUTED");
  expect(body.result.length).toBe(12);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(fasta);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("metric-length")).toContainText("12");
  await expect(page.getByTestId("metric-gc")).not.toContainText("N/A");
  await expect(page.getByTestId("invalid-input")).toHaveCount(0);
});

test("invalid DNA shows invalid-input not N/A metrics", async ({ page }) => {
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill("NNNN");
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("invalid-input")).toBeVisible();
  await expect(page.getByTestId("metric-grid")).toHaveCount(0);
  await expect(page.getByTestId("interpretation-panel")).toHaveCount(0);
});

test("AT-only DNA keeps GC 0 and does not turn NaN skew into 0 in the API", async ({ request }) => {
  const sequence = "AT".repeat(80);
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence, include_profiles: true, include_explanation: false },
  });
  const body = await api.json();
  expect(body.result.gc_content).toBe(0);
  expect(body.result.windowed_profiles[0].gc_skew.value_state).toBe("NAN");
  expect(body.result.windowed_profiles[0].gc_skew.value).toBeNull();
});

test("AT-only DNA UI shows Undefined skew not a fake zero", async ({ page }) => {
  const sequence = "AT".repeat(80);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-gc")).toContainText("0");
  await expect(page.getByTestId("metric-gc_skew")).toContainText("Undefined");
});

test("DNA ORF coordinates match API PREDICTED ORF", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: ORF_150, include_orfs: true, include_profiles: false, include_explanation: false, min_orf_length: 150 },
  });
  const body = await api.json();
  expect(body.result.orfs.length).toBeGreaterThan(0);
  const orf = body.result.orfs[0];
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("dna-orfs")).toContainText(String(orf.start));
  await expect(page.getByTestId("dna-orfs")).toContainText(String(orf.end));
  await expect(page.getByTestId("dna-orfs")).toContainText("PREDICTED");
});

test("DNA 50k length and GC match API", async ({ page, request }) => {
  const sequence = "A".repeat(50_000);
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence, include_orfs: false, include_profiles: false, include_explanation: false },
  });
  const body = await api.json();
  expect(body.result.length).toBe(50_000);
  await page.goto("/analysis/dna");
  await page.getByLabel("Include ORFs").uncheck();
  await page.getByLabel("Include windowed profiles").uncheck();
  await page.getByLabel("Include restriction sites").uncheck();
  await page.getByLabel("Include CpG islands").uncheck();
  await page.getByLabel("Include k-mers").uncheck();
  await page.getByLabel("Include SantaLucia Tm").uncheck();
  await fillControlled(page.locator("#dna-sequence"), sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-length")).toContainText("50000", { timeout: 60_000 });
  await expect(page.getByTestId("metric-gc")).toContainText(String(body.result.gc_content));
});
