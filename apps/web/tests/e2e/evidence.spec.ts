import { mkdirSync } from "node:fs";
import { resolve } from "node:path";
import { expect, test } from "@playwright/test";

const OUT = resolve(__dirname, "../../../../MIGRATION_evidence/prompt4_nextjs");

const PAGES = [
  ["overview", "/overview"],
  ["dna", "/analysis/dna"],
  ["rna", "/analysis/rna"],
  ["protein", "/analysis/protein"],
  ["alignment", "/analysis/alignment"],
  ["motif", "/discovery/motif"],
  ["ncbi", "/data/ncbi"],
  ["blast", "/data/blast"],
  ["msa", "/evolution/msa"],
  ["phylogeny", "/evolution/phylogeny"],
  ["crispr", "/discovery/crispr"],
  ["variant", "/discovery/variant"],
  ["engines", "/engines"],
  ["references", "/data/references"],
  ["viewer", "/structure/viewer"],
  ["compare", "/structure/compare"],
] as const;

test("capture workstation evidence screenshots", async ({ page }) => {
  mkdirSync(OUT, { recursive: true });
  await page.setViewportSize({ width: 1440, height: 900 });
  for (const [name, path] of PAGES) {
    await page.goto(path);
    await expect(page.locator(".hs-shell")).toBeVisible();
    await page.screenshot({ path: resolve(OUT, `${name}.png`), fullPage: true });
  }
  await page.goto("/structure/viewer");
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await page.screenshot({ path: resolve(OUT, "3d-1crn.png"), fullPage: true });
  await page.getByTestId("classify-4un3").click();
  await expect(page.getByTestId("mapping-status")).toContainText("UNCERTAIN");
  await page.screenshot({ path: resolve(OUT, "4UN3.png"), fullPage: true });
  await page.goto("/overview");
  await page.getByTestId("command-open").click();
  await expect(page.getByTestId("command-input")).toBeVisible();
  await page.screenshot({ path: resolve(OUT, "command-palette.png") });
});
