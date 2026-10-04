import { expect, test } from "@playwright/test";

test("BLAST availability declares local is not nt/nr; local job is engine-honest", async ({ page }) => {
  await page.goto("/data/blast");
  await expect(page.getByTestId("blast-not-nt")).toContainText("true");
  await expect(page.getByTestId("blast-local-status")).toBeVisible();
  await page.locator("#blast-q").fill("ATGCATGCATGC");
  await page.getByRole("button", { name: "Submit job" }).click();
  await expect(page.getByTestId("scientific-error").or(page.getByTestId("blast-hits")).or(page.getByTestId("job-state"))).toBeVisible({
    timeout: 60_000,
  });
  const err = page.getByTestId("scientific-error");
  if (await err.isVisible()) {
    await expect(err).toHaveAttribute("data-error-code", /ENGINE_NOT_INSTALLED|RESOURCE_LIMIT|ENGINE_FAILED|INVALID_INPUT/);
  }
});
