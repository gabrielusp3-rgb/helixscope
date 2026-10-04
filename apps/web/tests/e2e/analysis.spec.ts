import { expect, test } from "@playwright/test";
import { API } from "./helpers";

test("RNA analyze renders API CAI", async ({ page, request }) => {
  const sequence = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA";
  const api = await request.post(`${API}/api/v1/rna/analyze`, {
    data: { sequence, fold: false, include_codon_metrics: true, include_explanation: true },
  });
  const body = await api.json();
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText(String(body.result.length));
});

test("RNA transcribe replaces DNA T with U from the API", async ({ page, request }) => {
  const sequence = "ATGCATGC";
  const api = await request.post(`${API}/api/v1/rna/transcribe`, { data: { sequence } });
  const body = await api.json();
  expect(body.result.rna).toBe("AUGCAUGC");
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(sequence);
  await page.getByRole("button", { name: "Transcribe DNA to RNA" }).click();
  await expect(page.locator("#rna-seq")).toHaveValue(body.result.rna);
});

test("RNA fold is PREDICTED or ENGINE_NOT_INSTALLED, never invented", async ({ page, request }) => {
  const sequence = "GGGGAAACCCC";
  const api = await request.post(`${API}/api/v1/rna/analyze`, {
    data: { sequence, fold: true, include_codon_metrics: false, include_explanation: false },
  });
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(sequence);
  await page.getByLabel("Fold (ViennaRNA if installed)").check();
  await page.getByRole("button", { name: "Analyze" }).click();
  if (api.status() === 503) {
    await expect(page.getByTestId("scientific-error")).toHaveAttribute("data-error-code", "ENGINE_NOT_INSTALLED");
  } else {
    const body = await api.json();
    await expect(page.getByTestId("fold-status")).toContainText(String(body.result.fold.status));
    if (body.result.fold.mfe_kcal_mol != null) {
      await expect(page.getByTestId("rna-fold")).toContainText(String(body.result.fold.mfe_kcal_mol));
    }
  }
});

test("protein MW pI GRAVY match API", async ({ page, request }) => {
  const sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ";
  const api = await request.post(`${API}/api/v1/protein/analyze`, {
    data: { sequence, include_explanation: true },
  });
  const body = await api.json();
  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText(String(body.result.length));
  await expect(page.getByTestId("metric-pi")).toContainText(String(body.result.isoelectric_point).slice(0, 4));
  await expect(page.getByTestId("metric-gravy")).toContainText(String(body.result.gravy).slice(0, 4));
});

test("alignment identity matches API", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/alignment/pairwise`, {
    data: { seq1: "ACGTACGTACGT", seq2: "ACGTACGGACGT", mode: "global", include_explanation: true },
  });
  const body = await api.json();
  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("aligned-strings")).toContainText(body.result.aligned_seq1);
});

test("local Smith-Waterman strings match API", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/alignment/pairwise`, {
    data: { seq1: "ACGTACGTACGT", seq2: "ACGTACGGACGT", mode: "local", include_explanation: false },
  });
  const body = await api.json();
  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.locator("#mode").selectOption("local");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("metric-method")).toContainText("Smith-Waterman");
  await expect(page.getByTestId("aligned-strings")).toContainText(body.result.aligned_seq1);
});
