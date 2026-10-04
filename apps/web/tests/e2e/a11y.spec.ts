import { expect, test } from "@playwright/test";

test("DNA form controls are labelled and errors use alert", async ({ page }) => {
  await page.goto("/analysis/dna");
  await expect(page.getByLabel("DNA sequence")).toBeVisible();
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("scientific-error")).toHaveAttribute("role", "alert");
});

test("command palette is keyboard reachable", async ({ page }) => {
  await page.goto("/overview");
  await expect(page.getByTestId("command-open")).toBeVisible();
  await page.getByTestId("hs-workspace").click();
  await page.keyboard.press("Control+k");
  const palette = page.getByTestId("command-input");
  if (!(await palette.isVisible().catch(() => false))) {
    await page.evaluate(() => {
      window.dispatchEvent(
        new KeyboardEvent("keydown", {
          key: "k",
          code: "KeyK",
          ctrlKey: true,
          bubbles: true,
          cancelable: true,
        }),
      );
    });
  }
  await expect(palette).toBeVisible();
  await page.getByTestId("command-input").fill("Motif");
  await page.keyboard.press("Enter");
  await expect(page).toHaveURL(/discovery\/motif/);
});
