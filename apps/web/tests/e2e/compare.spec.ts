import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { expect, test } from "@playwright/test";

test("compare 8HSK vs 8HSF shows distinct RMSD fields from API", async ({ page, request }) => {
  const payload = JSON.parse(
    readFileSync(resolve(__dirname, "../../../../tests/fixtures/rcsb_alignment_complete.json"), "utf8"),
  );
  const api = await request.post("http://127.0.0.1:8000/api/v1/compare/parse", { data: { payload } });
  const body = await api.json();
  await page.goto("/structure/compare");
  await page.locator("#cmp-json").fill(JSON.stringify(payload));
  await page.getByRole("button", { name: "Parse alignment" }).click();
  const hero = page.getByTestId("compare-hero");
  await expect(hero).toContainText("Block RMSD");
  await expect(hero).toContainText("Global RMSD");
  await expect(hero).toContainText("Aligned residue pairs");
  await expect(hero).toContainText("Coverage");
  await expect(hero).toContainText("TM-score");
  await expect(hero).toContainText(String(body.result.n_aligned_residue_pairs));
  await expect(hero).toContainText(String(body.result.rmsd_block0_angstrom));
  await expect(hero).toContainText(String(body.result.rmsd_global_angstrom));
  await expect(page.getByTestId("compare-overlay-empty")).toBeVisible();
});
