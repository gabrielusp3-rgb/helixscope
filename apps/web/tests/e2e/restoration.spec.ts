import { expect, test } from "@playwright/test";
import { API } from "./helpers";

const DNA = "ATGC".repeat(40);
const RNA = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA";
const PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQKASEDLKKH";
const CRISPR = "ACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGT";

test("DNA restored surface shows metrics, charts, tables, interpretation, provenance", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: {
      sequence: DNA,
      include_orfs: true,
      include_profiles: true,
      include_restriction: true,
      include_cpg: true,
      include_kmers: true,
      include_santalucia: true,
      include_explanation: true,
    },
  });
  expect(api.ok()).toBeTruthy();
  const body = await api.json();
  expect(body.result.dinucleotides.CG).toBeTruthy();
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(DNA);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-length")).toContainText(String(body.result.length));
  await expect(page.getByTestId("dna-composition-bar")).toBeVisible();
  await expect(page.getByTestId("dna-dinuc-heatmap")).toBeVisible();
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
  await expect(page.getByTestId("dna-at-chart")).toBeVisible();
  await expect(page.getByTestId("dna-gc-skew-chart")).toBeVisible();
  await expect(page.getByTestId("dna-at-skew-chart")).toBeVisible();
  await expect(page.getByTestId("dna-entropy-chart")).toBeVisible();
  await expect(page.getByTestId("dna-kmer-bar")).toBeVisible();
  await expect(page.getByTestId("dna-sequence-viewer")).toBeVisible();
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
});

test("RNA restored surface shows composition, dinucleotide heatmap, RSCU", async ({ page }) => {
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(RNA);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();
  await expect(page.getByTestId("rna-composition")).toBeVisible();
  await expect(page.getByTestId("rna-dinuc-heatmap")).toBeVisible();
  await expect(page.getByTestId("rna-codon-heatmap")).toBeVisible();
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
});

test("Protein restored surface shows composition, radar, hydropathy", async ({ page }) => {
  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(PROTEIN);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();
  await expect(page.getByTestId("protein-aa-bar")).toBeVisible();
  await expect(page.getByTestId("protein-aa-radar")).toBeVisible();
  await expect(page.getByTestId("protein-hydropathy")).toBeVisible();
  await expect(page.getByTestId("protein-charge")).toBeVisible();
});

test("CRISPR restored surface shows guides and ranking, not only the form", async ({ page }) => {
  await page.goto("/discovery/crispr");
  await page.locator("#crispr-dna").fill(CRISPR);
  await page.getByLabel("Sequence-local off-target scan (not genome-wide)").check();
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guides")).toBeVisible();
  await expect(page.getByTestId("crispr-guide-scatter")).toBeVisible();
  await expect(page.getByTestId("crispr-ot-hits")).toBeVisible();
  await expect(page.getByTestId("crispr-primers")).toBeVisible();
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
});

test("Alignment restored surface shows column map and aligned strings", async ({ page }) => {
  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("aligned-strings")).toBeVisible();
  await expect(page.getByTestId("alignment-column-map")).toBeVisible();
});
