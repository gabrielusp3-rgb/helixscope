import { expect, test } from "@playwright/test";
import { API, ORF_150, PREALIGNED } from "./helpers";

const RNA = "AUGGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUGCUUAA";
const PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQKASEDLKKH";
const CRISPR = "ACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGT";

test("DNA FASTA yields computed metrics, features, infographics, explanation, raw data, provenance", async ({
  page,
  request,
}) => {
  const fasta = `>rebuild_dna description 9606\n${"ATGC".repeat(40)}\n`;
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: {
      sequence: fasta,
      include_orfs: true,
      include_profiles: true,
      include_restriction: true,
      include_cpg: true,
      include_kmers: true,
      include_explanation: true,
    },
  });
  expect(api.ok()).toBeTruthy();
  const body = await api.json();
  expect(body.result.status).toBe("COMPUTED");
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(fasta);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-status")).toContainText("COMPUTED");
  await expect(page.getByTestId("metric-length")).toContainText(String(body.result.length));
  await expect(page.getByTestId("metric-gc")).toContainText(String(body.result.gc_content));
  await expect(page.getByTestId("anonymous-sequence-notice")).toBeVisible();
  await expect(page.getByTestId("analysis-nav")).toBeVisible();
  await expect(page.getByTestId("dna-gc-chart")).toBeVisible();
  await expect(page.getByTestId("dna-entropy-chart")).toBeVisible();
  await expect(page.getByTestId("dna-composition-bar")).toBeVisible();
  await expect(page.getByTestId("what-this-means")).toBeVisible();
  await expect(page.getByTestId("method-limitations")).toBeVisible();
  await expect(page.getByTestId("metric-explanations")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
  await expect(page.getByTestId("interpretation-panel")).toBeVisible();
});

test("DNA ORF fixture exposes PREDICTED ORFs and feature table", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/dna/analyze`, {
    data: { sequence: ORF_150, include_orfs: true, include_profiles: false, min_orf_length: 150 },
  });
  const body = await api.json();
  expect(body.result.orfs.length).toBeGreaterThan(0);
  await page.goto("/analysis/dna");
  await page.locator("#dna-sequence").fill(ORF_150);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("dna-orfs")).toContainText("PREDICTED");
  await expect(page.getByTestId("dna-orfs")).toContainText(String(body.result.orfs[0].start));
});

test("RNA composition, codon metrics, explanation and provenance", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/rna/analyze`, {
    data: { sequence: RNA, fold: false, include_codon_metrics: true, include_explanation: true },
  });
  const body = await api.json();
  await page.goto("/analysis/rna");
  await page.locator("#rna-seq").fill(RNA);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText(String(body.result.length));
  await expect(page.getByTestId("rna-composition")).toBeVisible();
  await expect(page.getByTestId("rna-codon-heatmap")).toBeVisible();
  await expect(page.getByTestId("what-this-means")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
});

test("Protein FASTA exposes MW pI GRAVY composition and raw data", async ({ page, request }) => {
  const fasta = `>prot_fixture\n${PROTEIN}\n`;
  const api = await request.post(`${API}/api/v1/protein/analyze`, {
    data: { sequence: fasta, include_explanation: true },
  });
  const body = await api.json();
  await page.goto("/analysis/protein");
  await page.locator("#prot-seq").fill(fasta);
  await page.getByRole("button", { name: "Analyze" }).click();
  await expect(page.getByTestId("metric-len")).toContainText(String(body.result.length));
  await expect(page.getByTestId("metric-pi")).toContainText(String(body.result.isoelectric_point).slice(0, 4));
  await expect(page.getByTestId("metric-gravy")).toContainText(String(body.result.gravy).slice(0, 4));
  await expect(page.getByTestId("metric-inst")).toBeVisible();
  await expect(page.getByTestId("protein-aa-bar")).toBeVisible();
  await expect(page.getByTestId("what-this-means")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
});

test("Needleman-Wunsch and Smith-Waterman alignments render strings and column map", async ({ page, request }) => {
  const seq1 = "ACGTACGTACGT";
  const seq2 = "ACGTACGGACGT";
  const nw = await request.post(`${API}/api/v1/alignment/pairwise`, {
    data: { seq1, seq2, mode: "global", include_explanation: true, include_dotplot: true },
  });
  const nwBody = await nw.json();
  await page.goto("/analysis/alignment");
  await page.locator("#seq1").fill(seq1);
  await page.locator("#seq2").fill(seq2);
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("aligned-strings")).toContainText(nwBody.result.aligned_seq1);
  await expect(page.getByTestId("alignment-column-map")).toBeVisible();
  await expect(page.getByTestId("what-this-means")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
  await page.getByText("Edit alignment input").click();
  await page.locator("#mode").selectOption("local");
  await page.getByRole("button", { name: "Run" }).click();
  await expect(page.getByTestId("metric-method")).toContainText("Smith-Waterman");
});

test("CRISPR PAM locus yields guides, map, local scope, explanation", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/crispr/guides`, {
    data: { sequence: CRISPR, cas_system: "SpCas9", include_explanation: true, run_off_target: true },
  });
  const body = await api.json();
  expect(body.result.guides.length).toBeGreaterThan(0);
  await page.goto("/discovery/crispr");
  await page.locator("#crispr-dna").fill(CRISPR);
  await page.getByLabel("Sequence-local off-target scan (not genome-wide)").check();
  await page.getByRole("button", { name: "Design guides" }).click();
  await expect(page.getByTestId("crispr-guides")).toContainText(body.result.guides[0].guide_sequence);
  await expect(page.getByTestId("crispr-guides")).toContainText(String(body.result.guides[0].pam_sequence));
  await expect(page.getByTestId("crispr-guide-scatter")).toBeVisible();
  await expect(page.getByTestId("crispr-scope")).toBeVisible();
  await expect(page.getByTestId("what-this-means")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
  await expect(page.getByTestId("provenance-panel")).toBeVisible();
});

test("NCBI page is a gateway with PubMed option and no silent transfer", async ({ page }) => {
  await page.goto("/data/ncbi");
  await expect(page.getByTestId("external-disclosure")).toContainText("NCBI Entrez");
  await expect(page.locator("#db")).toContainText("pubmed");
  await expect(page.getByTestId("ncbi-workspace-transfer")).toHaveCount(0);
  await expect(page.getByTestId("ncbi-to-dna")).toHaveCount(0);
});

test("Scientific engines list remote providers without secrets", async ({ page }) => {
  await page.goto("/engines");
  await expect(page.getByTestId("remote-providers")).toBeVisible();
  await expect(page.getByTestId("remote-providers-table")).toContainText("NCBI");
  await expect(page.getByTestId("remote-providers-table")).not.toContainText("NCBI_API_KEY");
});

test("MSA conservation chart and phylogeny tree are Core-backed", async ({ page, request }) => {
  const msa = await request.post(`${API}/api/v1/msa/prealigned`, { data: { fasta: PREALIGNED } });
  const msaBody = await msa.json();
  await page.goto("/evolution/msa");
  await page.locator("#msa-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Import prealigned FASTA" }).click();
  await expect(page.getByTestId("msa-hash")).toContainText(msaBody.result.alignment_hash);
  await expect(page.getByTestId("msa-conservation")).toBeVisible();
  await expect(page.getByTestId("raw-data-panel")).toBeVisible();
  await page.goto("/evolution/phylogeny");
  await page.locator("#tree-fasta").fill(PREALIGNED);
  await page.getByRole("button", { name: "Infer tree" }).click();
  await expect(page.getByTestId("tree-viewer")).toBeVisible();
  await expect(page.getByTestId("newick")).not.toHaveText("");
});

test("Variant identity keeps exact integer coordinates", async ({ page, request }) => {
  const api = await request.post(`${API}/api/v1/variants/identify`, {
    data: { text: "17 43093557 C G", assembly: "GRCh38.p14" },
  });
  const body = await api.json();
  await page.goto("/discovery/variant");
  await page.locator("#vtext").fill("17 43093557 C G");
  await page.getByRole("button", { name: "Identify" }).click();
  await expect(page.getByTestId("metric-pos")).toContainText(String(body.result.position_0based));
  await expect(page.getByTestId("variant-position")).toContainText(String(body.result.position_0based));
  await expect(page.getByTestId("variant-position")).not.toContainText("e+");
});
