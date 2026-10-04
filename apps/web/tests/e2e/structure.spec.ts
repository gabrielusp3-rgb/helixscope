import { expect, test } from "@playwright/test";

test("3D fixtures load a canvas and 4UN3-class mapping stays UNCERTAIN", async ({ page, request }) => {
  const mapped = await request.post("http://127.0.0.1:8000/api/v1/structures/mapping/classify", {
    data: {
      sequences_identical: true,
      polymer_index_mode: false,
      observed_coordinate_residues: 0,
      polymer_length: 0,
    },
  });
  const mapBody = await mapped.json();
  expect(mapBody.result.mapping_status).toBe("UNCERTAIN");
  await page.goto("/structure/viewer");
  await page.getByRole("button", { name: "Load 1CRN" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible({ timeout: 60_000 });
  await page.getByRole("button", { name: "Load 1BNA" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
  await page.getByRole("button", { name: "Load 1RNA" }).click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
  await page.getByTestId("classify-4un3").click();
  await expect(page.getByTestId("mapping-status")).toContainText("UNCERTAIN");
  await expect(page.getByTestId("hand-status")).toContainText("OFF");
  await page.getByTestId("load-mapped-1crn").click();
  await expect(page.getByTestId("structure-viewport")).toBeVisible();
});
