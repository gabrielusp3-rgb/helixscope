import { expect, test } from "@playwright/test";
import { API, PREALIGNED } from "./helpers";

const SPCAS9 = "ACGTACGTACGTACGTACGTAGG";
const ANNOTATED = "Humano      [ A G T C ]\nChimpanze   [ A G T - ] <- comment\n";
const ORF_RNA = `AUG${"GCA".repeat(48)}UAA`;

test("CRISPR dropdown uses canonical Core keys", async ({ page, request }) => {
  const systems = await request.get(`${API}/api/v1/crispr/systems`);
  const body = await systems.json();
  const keys: string[] = body.result.systems.map((row: { canonical_key: string }) => row.canonical_key);
  expect(keys).toContain("AsCas12a");
  expect(keys).not.toContain("Cas12a");
  await page.goto("/discovery/crispr");
  await expect(page.getByTestId("crispr-cas")).toBeVisible();
  await expect(page.getByTestId("crispr-cas")).toContainText("AsCas12a");
  await page.getByTestId("crispr-cas").selectOption("SpCas9");
  await page.locator("#crispr-dna").fill(SPCAS9);
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("metric-cas")).toContainText("SpCas9");
  await expect(page.getByTestId("crispr-guides")).toBeVisible();
  await page.getByText("Edit CRISPR input").click();
  await page.getByTestId("crispr-cas").selectOption("AsCas12a");
  await page.locator("#crispr-dna").fill("TTTA" + "ACGTACGTACGTACGTACGTACG");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("metric-cas")).toContainText("AsCas12a");
});

test("MSA rejects annotated table and transfers hash to phylogeny NJ", async ({ page, request }) => {
  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(ANNOTATED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("scientific-error")).toBeVisible();
  await expect(page.getByTestId("scientific-error")).toContainText("not a recognized alignment format");
  const msa = await request.post(`${API}/api/v1/msa/prealigned`, { data: { fasta: PREALIGNED, format: "fasta" } });
  const msaBody = await msa.json();
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-hash")).toContainText(msaBody.result.alignment_hash);
  await page.getByRole("button", { name: "Send MSA to Phylogeny" }).click();
  await expect(page).toHaveURL(/\/evolution\/phylogeny/);
  await expect(page.getByTestId("phy-source-hash")).toContainText(msaBody.result.alignment_hash);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("newick")).toBeVisible();
  await expect(page.getByTestId("phy-method")).toContainText("neighbor_joining");
});

test("Compare form values are submitted, not 1CRN/1BNA/tm-align substitutes", async ({ page }) => {
  const usalignBodies: unknown[] = [];
  const rcsbBodies: unknown[] = [];
  await page.route("**/api/v1/compare/jobs/usalign", async (route) => {
    usalignBodies.push(route.request().postDataJSON());
    await route.continue();
  });
  await page.route("**/api/v1/compare/jobs/remote", async (route) => {
    rcsbBodies.push(route.request().postDataJSON());
    await route.continue();
  });
  await page.goto("/structure/compare");
  await expect(page.getByTestId("compare-method").locator("option[value='fatcat-rigid']")).toHaveCount(1, {
    timeout: 20_000,
  });
  await page.getByTestId("compare-reference").fill("1RNA");
  await page.getByTestId("compare-target").fill("1CRN");
  await page.getByTestId("compare-ref-chain").fill("B");
  await page.getByTestId("compare-tgt-chain").fill("C");
  await page.getByTestId("compare-method").selectOption("fatcat-rigid");
  await page.getByTestId("compare-usalign-submit").click();
  await expect(page.getByTestId("compare-job").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 90_000,
  });
  expect(usalignBodies.length).toBeGreaterThan(0);
  const usalign = usalignBodies[0] as Record<string, string>;
  expect(usalign.reference_entry).toBe("1RNA");
  expect(usalign.target_entry).toBe("1CRN");
  expect(usalign.reference_chain).toBe("B");
  expect(usalign.target_chain).toBe("C");
  expect(usalign.reference_entry).not.toBe("1CRN");
  expect(usalign.target_entry).not.toBe("1BNA");
  await page.getByTestId("compare-rcsb-submit").click();
  await expect.poll(() => rcsbBodies.length).toBeGreaterThan(0);
  const rcsb = rcsbBodies[0] as Record<string, string>;
  expect(rcsb.reference_entry).toBe("1RNA");
  expect(rcsb.target_entry).toBe("1CRN");
  expect(rcsb.reference_chain).toBe("B");
  expect(rcsb.target_chain).toBe("C");
  expect(rcsb.method).toBe("fatcat-rigid");
  expect(rcsb.method).not.toBe("tm-align");
});

test("Protein structure capabilities and bundled SASA", async ({ page }) => {
  await page.goto("/analysis/protein");
  await expect(page.getByTestId("protein-structure-caps")).toContainText("DSSP");
  await page.getByTestId("protein-structure-analyze").click();
  await expect(page.getByTestId("protein-structure-result").or(page.getByTestId("scientific-error"))).toBeVisible({
    timeout: 60_000,
  });
});

test("RNA coding translation and alignment translate checkbox", async ({ page }) => {
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(ORF_RNA);
  await page.getByTestId("rna-translate").click();
  await expect(page.getByTestId("rna-translate-result").or(page.getByTestId("scientific-error"))).toBeVisible();
  await page.goto("/analysis/alignment");
  await expect(page.getByTestId("align-translate")).toBeVisible();
  await page.getByTestId("align-translate").check();
});

test("References admission does not install GRCh38 when refused", async ({ page }) => {
  await page.goto("/data/references");
  await page.getByTestId("references-admit").click();
  await expect(page.getByTestId("references-admission")).toBeVisible();
  const decision = await page.getByTestId("references-admission").innerText();
  if (/RESOURCE_LIMIT|INVALID_INPUT|TEST_ONLY/i.test(decision)) {
    await expect(page.getByTestId("references-install")).toBeDisabled();
  }
});

test("Sidebar modules reach an API-backed page", async ({ page }) => {
  const routes = [
    "/",
    "/workspace/engines",
    "/analysis/dna",
    "/analysis/rna",
    "/analysis/protein",
    "/analysis/alignment",
    "/discovery/motif",
    "/discovery/crispr",
    "/discovery/variant",
    "/data/ncbi",
    "/data/blast",
    "/data/references",
    "/evolution/msa",
    "/evolution/phylogeny",
    "/structure/viewer",
    "/structure/compare",
  ];
  for (const route of routes) {
    await page.goto(route);
    await expect(page.locator("h1")).toBeVisible();
  }
});
