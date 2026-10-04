import type { Locator, Page } from "@playwright/test";

export const API = process.env.NEXT_PUBLIC_HELIX_API_BASE_URL ?? "http://127.0.0.1:8000";

export const PREALIGNED = `>seq_a
ACGTACGTACGTACGTACGT
>seq_b
ACGTACGTACGTACGTTTTT
>seq_c
TTTTACGTACGTACGTACGT
>seq_d
ACGTACGTAAAAACGTACGT
`;

export const ORF_150 = `ATG${"GCA".repeat(48)}TAA`;

export async function fillControlled(locator: Locator, value: string): Promise<void> {
  await locator.evaluate((el, seq) => {
    const tag = el.tagName;
    const proto =
      tag === "TEXTAREA" ? window.HTMLTextAreaElement.prototype : window.HTMLInputElement.prototype;
    const desc = Object.getOwnPropertyDescriptor(proto, "value");
    desc?.set?.call(el, seq);
    el.dispatchEvent(new Event("input", { bubbles: true }));
  }, value);
}

export async function waitScience(page: Page): Promise<void> {
  await page.getByTestId("metric-grid").or(page.getByTestId("scientific-error")).waitFor({ timeout: 60_000 });
}
