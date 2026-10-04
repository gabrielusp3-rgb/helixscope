import { expect, test, type Page } from "@playwright/test";
import { API } from "./helpers";

const INSULIN = `>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
`;

async function assertNoTopbarOverlap(page: Page) {
  const topbar = page.getByTestId("hs-topbar");
  const content = page.locator(".hs-page").first();
  await expect(topbar).toBeVisible();
  await expect(content).toBeVisible();
  const topBox = await topbar.boundingBox();
  const pageBox = await content.boundingBox();
  expect(topBox).toBeTruthy();
  expect(pageBox).toBeTruthy();
  if (!topBox || !pageBox) return;
  expect(topBox.y + topBox.height).toBeLessThanOrEqual(pageBox.y + 1);
}

test("API identity matches the frontend OpenAPI hash", async ({ page, request }) => {
  const identity = await request.get(`${API}/api/v1/system/identity`);
  expect(identity.ok()).toBeTruthy();
  const body = await identity.json();
  expect(body.openapi_sha256).toMatch(/^[a-f0-9]{64}$/);
  await page.goto("/overview");
  await expect(page.getByTestId("api-contract-mismatch")).toHaveCount(0);
  await expect(page.getByTestId("hs-topbar")).toContainText(String(body.product_version));
});

test("geometry: topbar does not overlap DNA page at 80/100/125 percent", async ({ page }) => {
  const viewports = [
    { width: 1366, height: 768 },
    { width: 1440, height: 900 },
    { width: 1920, height: 1080 },
  ];
  const zooms = [0.8, 1, 1.25];
  for (const viewport of viewports) {
    await page.setViewportSize(viewport);
    await page.goto("/analysis/dna");
    for (const zoom of zooms) {
      await page.evaluate((value) => {
        document.documentElement.style.zoom = String(value);
      }, zoom);
      await assertNoTopbarOverlap(page);
    }
  }
});

test("protein FASTA does not treat header punctuation as residues", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/protein/analyze`, {
    data: { sequence: INSULIN, include_explanation: true },
  });
  expect(api.ok()).toBeTruthy();
  const body = await api.json();
  expect(body.result.length).toBe(110);
  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(INSULIN);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText("110");
  await expect(page.getByTestId("protein-composition")).toBeVisible();
  await expect(page.getByTestId("scientific-error")).toHaveCount(0);
});

test("protein invalid residues are O and U only when those are the extras", async ({ request }) => {
  const api = await request.post(`${API}/api/v1/protein/analyze`, {
    data: { sequence: "MKTAYIOUAK", include_explanation: false },
  });
  expect(api.status()).toBe(400);
  const body = await api.json();
  expect(body.error.message).toMatch(/O, U/);
  expect(body.error.message).not.toContain(">");
  expect(body.error.message).not.toContain("|");
});

test("CRISPR 20-nt spacer is DNA but yields 0 guides with an explicit empty state", async ({ page }) => {
  await page.goto("/discovery/crispr");
  await page.locator("#crispr-dna").fill("ACGTACGTACGTACGTACGT");
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("metric-n")).toContainText("0");
  await expect(page.getByTestId("empty-state")).toContainText("20-nt spacer");
});

test("DNA send to CRISPR navigates and designs guides on a PAM locus", async ({ page, request }) => {
  const sequence = "ACGTACGTACGTACGTACGTAGG";
  const api = await request.post(`${API}/api/v1/crispr/guides`, {
    data: { sequence, cas_system: "SpCas9", include_explanation: true },
  });
  const body = await api.json();
  expect(body.result.n_guides).toBeGreaterThan(0);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(sequence);
  await page.getByRole("button", { name: "Send sequence to CRISPR" }).click();
  await expect(page).toHaveURL(/discovery\/crispr/);
  await expect(page.locator("#crispr-dna")).toHaveValue(sequence);
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guides")).toContainText(body.result.guides[0].guide_sequence);
});

test("DNA textarea value is the JSON sequence body", async ({ page }) => {
  const sequence = "ATGCATGCATGC";
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(sequence);
  const visible = await page.locator("#dna-sequence").inputValue();
  const pending = page.waitForRequest((req) => req.url().includes("/api/v1/dna/analyze") && req.method() === "POST");
  await page.getByRole("button", { name: "Analyze" }).click();
  const req = await pending;
  const payload = req.postDataJSON() as { sequence?: string };
  expect(payload.sequence).toBe(visible);
  expect(visible).toBe(sequence);
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("metric-length")).toContainText("12");
  await expect(page.getByTestId("dna-profiles-too-short")).toBeVisible();
});

test("DNA windowed graphs render for sequences at least 100 nt", async ({ page, request }) => {
  const sequence = "ATGC".repeat(40);
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence, include_profiles: true, include_explanation: false },
  });
  const body = await api.json();
  expect(body.result.windowed_profiles.length).toBeGreaterThan(0);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(sequence);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
  await expect(page.getByTestId("dna-at-chart")).toBeVisible();
  await expect(page.getByTestId("dna-gc-skew-chart")).toBeVisible();
  await expect(page.getByTestId("dna-at-skew-chart")).toBeVisible();
  await expect(page.getByTestId("dna-entropy-chart")).toBeVisible();
});

test("protein textarea value is the JSON sequence body", async ({ page }) => {
  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(INSULIN);
  const visible = await page.locator("#prot-seq").inputValue();
  const pending = page.waitForRequest((req) => req.url().includes("/api/v1/protein/analyze") && req.method() === "POST");
  await page.getByRole("button", { name: "Analyze" }).click();
  const req = await pending;
  const payload = req.postDataJSON() as { sequence?: string };
  expect(payload.sequence).toBe(visible);
  expect(payload.sequence).toContain(">");
  await expect(page.getByTestId("metric-len")).toContainText("110");
  await expect(page.getByTestId("scientific-error")).toHaveCount(0);
});

test("high-GC and low-GC DNA stay COMPUTED with finite GC", async ({ request }) => {
  const high = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: "GCGCGCGCGCGC", include_explanation: false, include_profiles: true },
  });
  const low = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: "ATATATATATAT", include_explanation: false, include_profiles: true },
  });
  const highBody = await high.json();
  const lowBody = await low.json();
  expect(highBody.result.status).toBe("COMPUTED");
  expect(lowBody.result.status).toBe("COMPUTED");
  expect(highBody.result.gc_content).toBe(100);
  expect(lowBody.result.gc_content).toBe(0);
});

test("geometry: protein and CRISPR pages stay below the topbar", async ({ page }) => {
  for (const route of ["/analysis/protein", "/discovery/crispr"]) {
    await page.setViewportSize({ width: 1366, height: 768 });
    await page.goto(route);
    await page.evaluate(() => {
      document.documentElement.style.zoom = "0.8";
    });
    await assertNoTopbarOverlap(page);
  }
});
