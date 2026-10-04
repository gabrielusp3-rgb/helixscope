import { expect, test } from "@playwright/test";
import { API } from "./helpers";

const VARIANT = "17 43093557 C G";
const ITEMS = [
  {
    field: "Variant A most_severe_consequence",
    value: "missense_variant",
    source: "Ensembl VEP",
    evidence_status: "RETRIEVED",
    retrieved_at_utc: "2026-08-29T00:00:00Z",
  },
  {
    field: "Variant B most_severe_consequence",
    value: "synonymous_variant",
    source: "NCBI ClinVar",
    evidence_status: "RETRIEVED",
    retrieved_at_utc: "2026-08-29T00:00:00Z",
  },
];

test("variant identify hash matches API and is not a diagnosis", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/variants/identify`, {
    data: { text: VARIANT, assembly: "GRCh38.p14" },
  });
  const body = await api.json();
  await page.goto("/discovery/variant");
  await page.locator("#vtext").fill(VARIANT);
  await page.getByRole("button", { name: "Identify" }).click();
  await expect(page.getByTestId("identity-hash")).toHaveText(body.result.identity_hash);
  await expect(page.getByTestId("metric-pos")).toContainText(String(body.result.position_0based));
  await expect(page.getByTestId("no-diagnosis")).toContainText("HelixScope does not say pathogenic");
});

test("variant explore with remotes off keeps effect NOT_COMPUTED", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/variants/explore`, {
    data: { text: VARIANT, assembly: "GRCh38.p14", enable_vep: false, enable_clinvar: false, enable_domains: false },
  });
  const body = await api.json();
  const identity = body.result.variant ?? body.result.layers?.variant_identity?.data ?? body.result;
  expect(identity.effect_status).toBe("NOT_COMPUTED");
  await page.goto("/discovery/variant");
  await page.locator("#vtext").fill(VARIANT);
  await page.getByRole("button", { name: "Explore", exact: true }).click();
  await expect(page.getByTestId("metric-effect")).toContainText("NOT_COMPUTED");
  await expect(page.getByTestId("no-diagnosis")).toBeVisible();
});

test("evidence pack keeps confidence null and lists conflicts", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/evidence/pack`, {
    data: { kind: "variants", items: ITEMS },
  });
  const body = await api.json();
  expect(body.result.confidence_score).toBeNull();
  await page.goto("/discovery/variant");
  await page.locator("#evidence-json").fill(JSON.stringify(ITEMS));
  await page.getByRole("button", { name: "Pack evidence" }).click();
  await expect(page.getByTestId("metric-confidence")).toContainText("N/A");
  await expect(page.getByTestId("evidence-conflicts")).toBeVisible();
});
