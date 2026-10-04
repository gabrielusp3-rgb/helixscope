import { expect, test } from "@playwright/test";

test("CRISPR starts empty then fixture guides match API null advanced models", async ({ page, request }) => {
  await page.goto("/discovery/crispr");
  await expect(page.locator("#crispr-dna")).toHaveValue("");
  await expect(page.getByTestId("crispr-name")).toHaveValue("");
  const sequence = "ACGTACGTACGTACGTACGTAGG";
  const api = await request.post("http://127.0.0.1:8000/api/v1/crispr/guides", {
    data: { sequence, cas_system: "SpCas9", include_explanation: true },
  });
  const body = await api.json();
  await page.locator("#crispr-dna").fill(sequence);
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guides")).toContainText(body.result.guides[0].guide_sequence);
  await expect(page.getByTestId("crispr-scope")).toContainText("TEST REFERENCE");
  expect(body.result.guides[0].ruleset2_score ?? body.result.guides[0].deephf_score ?? null).toBeNull();
});
