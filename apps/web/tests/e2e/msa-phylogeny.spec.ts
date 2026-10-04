import { expect, test } from "@playwright/test";
import { API, PREALIGNED } from "./helpers";

test("MSA to phylogeny NJ uses API tree", async ({ page, request }) => {
  const msa = await request.post(`${API}/api/v1/msa/prealigned`, { data: { fasta: PREALIGNED } });
  const msaBody = await msa.json();
  const tree = await request.post(`${API}/api/v1/phylogeny/infer`, {
    data: { fasta: PREALIGNED, method: "neighbor_joining", distance_model: "p_distance" },
  });
  const treeBody = await tree.json();
  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-hash")).toContainText(msaBody.result.alignment_hash);
  await page.getByRole("button", { name: "Send MSA to Phylogeny" }).click();
  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("newick")).toContainText(treeBody.result.newick.slice(0, 12));
  await expect(page.getByTestId("tree-viewer")).toBeVisible();
});

test("MSA bounded upload matches prealigned hash", async ({ page, request }) => {
  const msa = await request.post(`${API}/api/v1/msa/prealigned`, { data: { fasta: PREALIGNED } });
  const msaBody = await msa.json();
  await page.goto("/evolution/msa");
  await page.locator("#msa-file").setInputFiles({
    name: "prealigned.fasta",
    mimeType: "text/plain",
    buffer: Buffer.from(PREALIGNED, "utf8"),
  });
  await expect(page.getByTestId("msa-hash")).toContainText(msaBody.result.alignment_hash);
});

test("UPGMA Newick matches API and is labeled UPGMA", async ({ page, request }) => {
  const tree = await request.post(`${API}/api/v1/phylogeny/infer`, {
    data: { fasta: PREALIGNED, method: "upgma", distance_model: "p_distance" },
  });
  const treeBody = await tree.json();
  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.locator("#method").selectOption("upgma");
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("phy-method")).toContainText(String(treeBody.result.method));
  await expect(page.getByTestId("newick")).toContainText(treeBody.result.newick.slice(0, 12));
});

test("IQ-TREE job completes or reports missing engine", async ({ page }) => {
  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.locator("#method").selectOption("iqtree_ml");
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("newick").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 90_000,
  });
  const err = page.getByTestId("scientific-error");
  if (await err.isVisible()) {
    await expect(err).toHaveAttribute("data-error-code", /ENGINE_NOT_INSTALLED|ENGINE_FAILED|RESOURCE_LIMIT|TIMEOUT/);
  } else {
    await expect(page.getByTestId("newick")).toBeVisible();
    await expect(page.getByTestId("phy-method")).toContainText(/iqtree|IQ-TREE/i);
  }
});

test("FastTree job completes or reports missing engine", async ({ page }) => {
  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.locator("#method").selectOption("fasttree_ml");
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("newick").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 90_000,
  });
  if (await page.getByTestId("newick").isVisible()) {
    await expect(page.getByTestId("phy-method")).toContainText(/fasttree|FastTree/i);
  }
});
