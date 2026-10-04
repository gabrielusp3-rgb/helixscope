import { expect, test } from "@playwright/test";
import { API, ORF_150, PREALIGNED } from "./helpers";

test("DNA workstation shows computed metrics, dark chart, ORF track and provenance", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: ORF_150, include_orfs: true, include_profiles: true, include_explanation: true, min_orf_length: 150 },
  });
  const body = await api.json();
  expect(body.result.status).toBe("COMPUTED");
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("metric-gc")).toContainText(String(body.result.gc_content));
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
  await expect(page.getByTestId("dna-dinuc-heatmap")).toBeVisible();
  await expect(page.getByTestId("genomic-track")).toBeVisible();
  await expect(page.getByTestId("dna-orfs")).toContainText("PREDICTED");
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
});

test("CRISPR Cas registry drives SpCas9 and AsCas12a locus inspector", async ({ page, request }) => {
  const systems = await request.get(`${API}/api/v1/crispr/systems`);
  const registry = await systems.json();
  const keys: string[] = registry.result.systems.map((row: { canonical_key: string }) => row.canonical_key);
  expect(keys).toContain("AsCas12a");
  expect(keys).toContain("Cas13");
  await page.goto("/discovery/crispr");
  await expect(page.getByTestId("cas-system-card")).toBeVisible();
  await page.getByTestId("crispr-cas").selectOption("SpCas9");
  await expect(page.getByTestId("cas-canonical")).toContainText("SpCas9");
  await page.locator("#crispr-dna").fill("ACGTACGTACGTACGTACGTAGG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("metric-cas")).toContainText("SpCas9");
  await expect(page.getByTestId("crispr-locus")).toBeVisible();
  await expect(page.getByTestId("crispr-guide-inspector")).toBeVisible();
  await expect(page.getByTestId("crispr-scope")).toContainText("SEQUENCE LOCAL");
  await page.getByText("Edit CRISPR input").click();
  await page.getByTestId("crispr-cas").selectOption("AsCas12a");
  await expect(page.getByTestId("cas-canonical")).toContainText("AsCas12a");
  await expect(page.getByTestId("cas-pam")).toContainText("TTTV");
  await expect(page.getByTestId("cas-target-molecule")).toContainText("DNA");
  await page.locator("#crispr-dna").fill("TTTA" + "ACGTACGTACGTACGTACGTACG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("metric-cas")).toContainText("AsCas12a");
  await expect(page.getByTestId("cas-mode")).not.toContainText("identical-to-SpCas9");
});

test("Cas13 card stays RNA-aware before design", async ({ page }) => {
  await page.goto("/discovery/crispr");
  await expect(page.getByTestId("crispr-cas")).toContainText("Cas13");
  await page.getByTestId("crispr-cas").selectOption("Cas13");
  await expect(page.getByTestId("cas-target-molecule")).toContainText("RNA");
  await expect(page.getByTestId("cas-rna-note")).toBeVisible();
  await expect(page.getByLabel("Target RNA")).toBeVisible();
});

test("MSA viewer and phylogeny NJ remain typed", async ({ page, request }) => {
  const msa = await request.post(`${API}/api/v1/msa/prealigned`, { data: { fasta: PREALIGNED, format: "fasta" } });
  const msaBody = await msa.json();
  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-viewer")).toBeVisible();
  await expect(page.getByTestId("msa-hash")).toContainText(msaBody.result.alignment_hash);
  await page.getByRole("button", { name: "Send MSA to Phylogeny" }).click();
  await expect(page).toHaveURL(/evolution\/phylogeny/);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("tree-viewer")).toBeVisible();
  await expect(page.getByTestId("phy-method")).toContainText("neighbor_joining");
  await expect(page.getByTestId("taxonomy-not-phylogeny")).toContainText("TAXONOMIC ANNOTATION");
});

test("3D viewer loads bundled 1CRN with experimental banner", async ({ page }) => {
  await page.goto("/structure/viewer");
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible({ timeout: 60_000 });
  await expect(page.getByTestId("structure-evidence-banner")).toContainText("EXPERIMENTAL");
  await expect(page.getByTestId("structure-evidence-banner")).toContainText("1CRN");
  await page.getByTestId("classify-4un3").click();
  await expect(page.getByTestId("mapping-status")).toContainText("UNCERTAIN");
});

test("command palette is content-sized and presentation toggle exists", async ({ page }) => {
  await page.goto("/overview");
  await expect(page.getByTestId("presentation-mode")).toHaveValue("beginner");
  await page.getByTestId("presentation-mode").selectOption("expert");
  await expect(page.getByTestId("presentation-mode")).toHaveValue("expert");
  await page.getByTestId("command-open").click();
  const dialog = page.locator(".cmd-dialog");
  await expect(dialog).toBeVisible();
  const box = await dialog.boundingBox();
  expect(box).toBeTruthy();
  expect(box!.height).toBeLessThan(700);
  await page.getByTestId("command-input").fill("Phylogeny");
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/evolution\/phylogeny/);
});

test("protein bundled SASA panel remains reachable", async ({ page }) => {
  await page.goto("/analysis/protein");
  await page.getByTestId("protein-structure-search").click();
  await expect(page.getByTestId("protein-structure-result").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 60_000,
  });
});
