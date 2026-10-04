import { expect, test, type Page } from "@playwright/test";
import { ORF_150, PREALIGNED } from "./helpers";

const RNA = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA";
const PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQKASEDLKKH";

const ROUTES = [
  "/overview",
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
  "/engines",
];

type Capture = {
  consoleErrors: string[];
  pageErrors: string[];
  failedRequests: string[];
  httpErrors: string[];
};

function attachCapture(page: Page): Capture {
  const capture: Capture = { consoleErrors: [], pageErrors: [], failedRequests: [], httpErrors: [] };
  page.on("console", (msg) => {
    if (msg.type() === "error") capture.consoleErrors.push(msg.text());
  });
  page.on("pageerror", (err) => {
    capture.pageErrors.push(String(err));
  });
  page.on("requestfailed", (req) => {
    const url = req.url();
    if (url.includes("vercel.live") || url.includes("va.vercel-scripts.com")) return;
    capture.failedRequests.push(`${req.failure()?.errorText || "failed"} ${url}`);
  });
  page.on("response", (res) => {
    const status = res.status();
    const url = res.url();
    if (status < 400) return;
    if (url.includes("vercel.live") || url.includes("va.vercel-scripts.com")) return;
    capture.httpErrors.push(`${status} ${url}`);
  });
  return capture;
}

function assertNoFatal(capture: Capture): void {
  const fatalConsole = capture.consoleErrors.filter(
    (text) =>
      !text.includes("Report Only") &&
      !text.includes("Content-Security-Policy") &&
      !text.includes("vercel.live") &&
      !text.includes("Download the React DevTools"),
  );
  expect(capture.pageErrors, `pageerrors: ${capture.pageErrors.join(" | ")}`).toEqual([]);
  expect(fatalConsole, `console: ${fatalConsole.join(" | ")}`).toEqual([]);
  expect(capture.httpErrors, `http: ${capture.httpErrors.join(" | ")}`).toEqual([]);
}

async function waitContract(page: Page): Promise<void> {
  await expect(page.getByTestId("hs-topbar")).toHaveAttribute("data-contract-status", "ok", { timeout: 45_000 });
}

test("public routes load without fatal console or HTTP errors", async ({ page }) => {
  const capture = attachCapture(page);
  for (const route of ROUTES) {
    const response = await page.goto(route, { waitUntil: "domcontentloaded" });
    expect(response?.ok(), route).toBeTruthy();
    await expect(page.locator("h1").first()).toBeVisible();
  }
  assertNoFatal(capture);
});

test("public DNA RNA protein MSA 3D browser smoke", async ({ page }) => {
  const capture = attachCapture(page);

  await page.goto("/analysis/dna");
  await waitContract(page);
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");

  await page.goto("/analysis/rna");
  await waitContract(page);
  await page.locator("#rna-seq").fill(RNA);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toBeVisible();

  await page.goto("/analysis/protein");
  await waitContract(page);
  await page.locator("#prot-seq").fill(PROTEIN);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-grid").first()).toBeVisible();

  await page.goto("/evolution/msa");
  await waitContract(page);
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-viewer").or(page.getByTestId("metric-grid")).or(page.getByTestId("scientific-error"))).toBeVisible();

  await page.goto("/evolution/phylogeny");
  await waitContract(page);
  await expect(page.locator("h1")).toContainText("Phylogeny");

  await page.goto("/structure/viewer");
  await waitContract(page);
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
  const molstar = page.getByTestId("molstar-viewport");
  const molstarError = page.getByTestId("molstar-error");
  await expect(molstar.or(molstarError)).toBeVisible({ timeout: 60_000 });

  assertNoFatal(capture);
});
