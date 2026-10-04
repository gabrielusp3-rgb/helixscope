import { expect, test } from "@playwright/test";
import { API } from "./helpers";

test("NCBI fetch without accession surfaces INVALID_INPUT without contacting a live record", async ({ page }) => {
  await page.goto("/data/ncbi");
  await expect(page.getByTestId("external-disclosure")).toContainText("NCBI Entrez");
  await page.getByRole("button", { name: "Fetch record" }).click();
  await expect(page.getByTestId("scientific-error")).toBeVisible();
  await expect(page.getByTestId("scientific-error")).toHaveAttribute("data-error-code", "INVALID_INPUT");
});

test("NCBI search rejects URL-like terms locally", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/ncbi/search`, {
    data: { term: "http://example.invalid", email: "helixscope@example.invalid", db: "nucleotide" },
  });
  expect(api.ok()).toBeFalsy();
  await page.goto("/data/ncbi");
  await page.locator("#email").fill("helixscope@example.invalid");
  await page.locator("#term").fill("http://example.invalid");
  await page.getByRole("button", { name: "Search" }).click();
  await expect(page.getByTestId("scientific-error")).toBeVisible();
});
