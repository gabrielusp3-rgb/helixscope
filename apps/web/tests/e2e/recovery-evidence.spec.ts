import { mkdirSync, readFileSync } from "node:fs";
import { resolve } from "node:path";
import { expect, test, type Page } from "@playwright/test";
import { PREALIGNED } from "./helpers";

const OUT = resolve(process.cwd(), "../../FUNCTIONAL_RECOVERY_evidence");

const INSULIN = `>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
`;

async function shot(page: Page, name: string) {
  mkdirSync(OUT, { recursive: true });
  const file = resolve(OUT, `${name}.png`);
  await page.screenshot({ path: file, fullPage: false });
}

async function analyzeDna(page: Page) {
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill("ATGC".repeat(40));
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
}

test("capture functional recovery evidence after real analyses", async ({ page }) => {
  mkdirSync(OUT, { recursive: true });

  const dnaViewports: Array<{ name: string; width: number; height: number; zoom: number }> = [
    { name: "dna-1366x768-80", width: 1366, height: 768, zoom: 0.8 },
    { name: "dna-1366x768-100", width: 1366, height: 768, zoom: 1 },
    { name: "dna-1440x900-100", width: 1440, height: 900, zoom: 1 },
    { name: "dna-1920x1080-100", width: 1920, height: 1080, zoom: 1 },
  ];
  for (const viewport of dnaViewports) {
    await page.setViewportSize({ width: viewport.width, height: viewport.height });
    await page.evaluate((value) => {
      document.documentElement.style.zoom = String(value);
    }, viewport.zoom);
    await analyzeDna(page);
    await expect(page.getByTestId("hs-topbar")).toBeVisible();
    await expect(page.getByRole("heading", { name: "DNA", exact: true })).toBeVisible();
    await shot(page, viewport.name);
  }

  await page.setViewportSize({ width: 1920, height: 1080 });
  await page.evaluate(() => {
    document.documentElement.style.zoom = "1";
  });

  await analyzeDna(page);
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
  await shot(page, "dna-graphs");
  await page.getByRole("button", { name: "Load illustrative B-DNA helix" }).click();
  await expect(page.getByTestId("dna-3d-kind")).toContainText("ILLUSTRATIVE");
  await shot(page, "dna-3d");

  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill("AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA");
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();
  await shot(page, "rna-result");

  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(INSULIN);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText("110");
  await shot(page, "protein-result");

  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill("ACGTACGTACGT");
  await page.locator("#seq2").fill("ACGTACGGACGT");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("aligned-strings")).toBeVisible();
  await shot(page, "alignment");

  await page.goto("/discovery/crispr");
  await page.locator("#crispr-dna").fill("ACGTACGTACGTACGTACGTAGG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guides")).toBeVisible();
  await shot(page, "crispr");

  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-hash")).toBeVisible();
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
  await page.getByRole("button", { name: "Load 1RNA" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
  await shot(page, "3d-1rna");
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
  await shot(page, "3d-1crn");
  await page.getByTestId("classify-4un3").click();
  await expect(page.getByTestId("mapping-status")).toContainText("UNCERTAIN");
  await shot(page, "3d-4un3-uncertain");

  const payload = JSON.parse(
    readFileSync(resolve(__dirname, "../../../../tests/fixtures/rcsb_alignment_complete.json"), "utf8"),
  );
  await page.goto("/structure/compare");
  await page.locator("#cmp-json").fill(JSON.stringify(payload));
  await page.getByRole("button", { name: "Parse alignment" }).click();
  await expect(page.getByTestId("compare-hero")).toBeVisible();
  await shot(page, "compare");
});
