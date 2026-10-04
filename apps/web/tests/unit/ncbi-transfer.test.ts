import { describe, expect, it } from "vitest";
import {
  ncbiIsPubmed,
  ncbiRetrievedSequence,
  suggestedNcbiTarget,
} from "@/lib/workspace/ncbiTransfer";

describe("NCBI workspace transfer helpers", () => {
  it("does not invent a sequence from PubMed metadata", () => {
    const result = {
      record: { record_kind: "pubmed", title: "A paper", sequence: "" },
    };
    expect(ncbiIsPubmed(result, "pubmed")).toBe(true);
    expect(ncbiRetrievedSequence(result)).toBe("");
    expect(suggestedNcbiTarget(result, "pubmed")).toBeNull();
  });

  it("routes protein database records to protein analysis", () => {
    const result = {
      record: { accession: "P01308", sequence: "MALWMRLLPLLALLALWGPDPAAA", molecule: "protein" },
    };
    expect(suggestedNcbiTarget(result, "protein")).toBe("protein");
    expect(ncbiRetrievedSequence(result)).toBe("MALWMRLLPLLALLALWGPDPAAA");
  });

  it("defaults nucleotide DNA-like sequence to DNA, not a guessed gene", () => {
    const result = {
      record: { accession: "NM_000207", sequence: "ATGCATGC", molecule: "mRNA" },
    };
    expect(suggestedNcbiTarget(result, "nucleotide")).toBe("dna");
  });
});
