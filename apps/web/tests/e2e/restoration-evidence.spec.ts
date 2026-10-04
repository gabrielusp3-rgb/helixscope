import { mkdirSync } from "node:fs";
import { resolve } from "node:path";
import { expect, test, type Page } from "@playwright/test";
import { PREALIGNED } from "./helpers";

const OUT = resolve(process.cwd(), "../../TOTAL_RESTORATION_evidence/next_after");

const INSULIN = `>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
`;

async function shot(page: Page, name: string) {
  mkdirSync(OUT, { recursive: true });
  await page.screenshot({ path: resolve(OUT, `${name}.png`), fullPage: true });
}

test("capture restored Next.js scientific surfaces after real analyses", async ({ page }) => {
  mkdirSync(OUT, { recursive: true });
  await page.setViewportSize({ width: 1920, height: 1080 });

  await page.goto("/overview");
  await shot(page, "overview");

  await page.goto("/engines");
  await shot(page, "engines");

  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill("ATGC".repeat(40));
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("dna-dinuc-heatmap")).toBeVisible();
  await shot(page, "dna");
  await page.getByRole("button", { name: "Load illustrative B-DNA helix" }).click();
  await expect(page.getByTestId("dna-3d-kind")).toBeVisible();
  await shot(page, "dna-3d");

  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill("AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA");
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("rna-dinuc-heatmap")).toBeVisible();
  await shot(page, "rna");

  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(INSULIN);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("protein-aa-radar")).toBeVisible();
  await shot(page, "protein");

  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("alignment-column-map")).toBeVisible();
  await shot(page, "alignment");

  await page.goto("/discovery/motif");
  await page.locator("#motif-seq").fill("ACGTACGTACGT");
  await page.locator("#pattern").fill("ACGT");
  await page.getByRole("button", { name: "Search" }).click();
  await expect(page.getByTestId("motif-hits")).toBeVisible();
  await shot(page, "motif");

  await page.goto("/discovery/crispr");
  await page.locator("#crispr-dna").fill("ACGTACGTACGTACGTACGTAGG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guide-scatter")).toBeVisible();
  await shot(page, "crispr");

  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-viewer")).toBeVisible();
  await shot(page, "msa");

  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("tree-viewer")).toBeVisible();
  await shot(page, "phylogeny");

  await page.goto("/structure/viewer");
  await page.getByRole("button", { name: "Load 1BNA" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible({ timeout: 60_000 });
  await shot(page, "3d-1bna");
});
