import { mkdirSync } from "node:fs";
import { resolve } from "node:path";
import { expect, test, type Page } from "@playwright/test";
import { ORF_150, PREALIGNED } from "./helpers";

const ROOT = resolve(process.cwd(), "../../PROMPT2_evidence");

async function shot(page: Page, folder: string, name: string) {
  const dir = resolve(ROOT, folder);
  mkdirSync(dir, { recursive: true });
  await page.screenshot({ path: resolve(dir, `${name}.png`), fullPage: true });
}

test("evidence: overview engines DNA RNA protein alignment", async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/overview");
  await expect(page.getByRole("heading", { name: "HelixScope" })).toBeVisible();
  await shot(page, "overview", "overview");
  await page.goto("/engines");
  await shot(page, "overview", "engines");

  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("genomic-track")).toBeVisible();
  await shot(page, "dna", "overview");
  await page.getByTestId("dna-dinuc-heatmap").scrollIntoViewIfNeeded();
  await shot(page, "dna", "dinucleotide");
  await page.getByTestId("dna-gc-chart").scrollIntoViewIfNeeded();
  await shot(page, "dna", "profiles");
  await page.getByTestId("genomic-track").scrollIntoViewIfNeeded();
  await shot(page, "dna", "orf-track");
  await page.getByTestId("interpretation-panel").scrollIntoViewIfNeeded();
  await shot(page, "dna", "interpretation");
  await page.getByTestId("provenance-panel").scrollIntoViewIfNeeded();
  await shot(page, "dna", "provenance");
  await shot(page, "dna", "composition");

  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill("AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA");
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();
  await shot(page, "rna", "overview");

  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN");
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();
  await shot(page, "protein", "physicochemical");
  await page.getByTestId("protein-structure-search").click();
  await expect(page.getByTestId("protein-structure-result").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 60_000,
  });
  await shot(page, "protein", "structure-search");

  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("aligned-strings")).toBeVisible();
  await shot(page, "alignment", "pairwise");
});

test("evidence: CRISPR SpCas9 and AsCas12a", async ({ page }) => {
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/discovery/crispr");
  await expect(page.getByTestId("cas-system-card")).toBeVisible();
  await page.getByTestId("crispr-cas").selectOption("SpCas9");
  await page.locator("#crispr-dna").fill("ACGTACGTACGTACGTACGTAGG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-locus")).toBeVisible();
  await shot(page, "crispr", "spcas9-locus");
  await page.getByTestId("crispr-guide-inspector").scrollIntoViewIfNeeded();
  await shot(page, "crispr", "guide-inspector");
  await page.getByText("Edit CRISPR input").click();
  await page.getByTestId("crispr-cas").selectOption("AsCas12a");
  await expect(page.getByTestId("cas-pam")).toContainText("TTTV");
  await shot(page, "crispr", "ascas12a-card");
});

test("evidence: variant NCBI BLAST MSA phylogeny compare", async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/discovery/variant");
  await shot(page, "variant", "workspace");
  await page.goto("/data/ncbi");
  await shot(page, "ncbi", "workspace");
  await page.goto("/data/blast");
  await expect(page.getByTestId("blast-backend")).toBeVisible();
  await shot(page, "blast", "workspace");

  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-viewer")).toBeVisible();
  await shot(page, "msa", "viewer");
  await page.getByRole("button", { name: "Send MSA to Phylogeny" }).click();
  await expect(page).toHaveURL(/evolution\/phylogeny/);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("tree-viewer")).toBeVisible();
  await shot(page, "phylogeny", "nj");
  await page.locator("details.hs-input-compact").evaluate((el) => {
    (el as HTMLDetailsElement).open = true;
  });
  await page.locator("#method").selectOption("upgma");
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("phy-method")).toContainText("upgma");
  await shot(page, "phylogeny", "upgma");

  await page.goto("/structure/compare");
  await shot(page, "compare", "workspace");
});

test("evidence: experimental 3D banners", async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/structure/viewer");
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await expect(page.getByTestId("structure-evidence-banner")).toContainText("EXPERIMENTAL", { timeout: 60_000 });
  await shot(page, "structure", "1crn");
  await page.getByRole("button", { name: "Load 1BNA" }).click();
  await expect(page.getByTestId("structure-evidence-banner")).toContainText("1BNA", { timeout: 60_000 });
  await shot(page, "structure", "1bna");
  await page.getByRole("button", { name: "Load 1RNA" }).click();
  await expect(page.getByTestId("structure-evidence-banner")).toContainText("1RNA", { timeout: 60_000 });
  await shot(page, "structure", "1rna");
  await page.getByTestId("classify-4un3").click();
  await expect(page.getByTestId("mapping-status")).toContainText("UNCERTAIN");
  await shot(page, "structure", "4un3-uncertain");
});

test("evidence: illustrative B-DNA and k-mers", async ({ page }) => {
  test.setTimeout(120_000);
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  const kmers = page.getByTestId("dna-kmers");
  if (await kmers.count()) {
    await kmers.scrollIntoViewIfNeeded();
    await shot(page, "dna", "kmers");
  }
  await page.getByRole("button", { name: "Load illustrative B-DNA helix" }).click();
  await expect(page.getByTestId("dna-3d-kind")).toContainText("ILLUSTRATIVE");
  await shot(page, "structure", "illustrative-dna");
});

test("evidence: responsive overview and DNA", async ({ page }) => {
  await page.setViewportSize({ width: 1366, height: 768 });
  await page.goto("/overview");
  await shot(page, "responsive", "1366x768-100");
  await page.goto("/analysis/dna");
  await shot(page, "responsive", "dna-1366x768");
  await page.setViewportSize({ width: 1440, height: 900 });
  await page.goto("/overview");
  await shot(page, "responsive", "1440x900-100");
  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.goto("/overview");
  await shot(page, "responsive", "1920x1080-100");
});
