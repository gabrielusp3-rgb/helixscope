export type NavItem = {
  href: string;
  label: string;
  group: string;
  keywords: string[];
};

export const NAV_GROUPS = [
  "Workspace",
  "Analysis",
  "Discovery",
  "Data",
  "Evolution",
  "Structure",
] as const;

export const NAV_ITEMS: NavItem[] = [
  { href: "/overview", label: "Overview", group: "Workspace", keywords: ["home", "start"] },
  { href: "/engines", label: "Scientific Engines", group: "Workspace", keywords: ["tools", "status"] },
  { href: "/analysis/dna", label: "DNA", group: "Analysis", keywords: ["gc", "orf"] },
  { href: "/analysis/rna", label: "RNA", group: "Analysis", keywords: ["fold", "cai", "rscu"] },
  { href: "/analysis/protein", label: "Protein", group: "Analysis", keywords: ["gravy", "pi"] },
  { href: "/analysis/alignment", label: "Alignment", group: "Analysis", keywords: ["needleman", "smith"] },
  { href: "/discovery/motif", label: "Motif Search", group: "Discovery", keywords: ["iupac"] },
  { href: "/discovery/crispr", label: "CRISPR Design", group: "Discovery", keywords: ["cas9", "guide"] },
  { href: "/discovery/variant", label: "Variant Explorer", group: "Discovery", keywords: ["vep", "clinvar"] },
  { href: "/data/ncbi", label: "NCBI Fetch", group: "Data", keywords: ["entrez", "accession"] },
  { href: "/data/blast", label: "NCBI BLAST", group: "Data", keywords: ["blast+", "remote"] },
  { href: "/data/references", label: "References", group: "Data", keywords: ["grch38", "t2t"] },
  { href: "/evolution/msa", label: "MSA", group: "Evolution", keywords: ["alignment", "clustal"] },
  { href: "/evolution/phylogeny", label: "Phylogeny", group: "Evolution", keywords: ["tree", "iqtree"] },
  { href: "/structure/viewer", label: "3D Viewer", group: "Structure", keywords: ["pdb", "1crn"] },
  { href: "/structure/compare", label: "Structure Compare", group: "Structure", keywords: ["rmsd", "tm-score"] },
];

export function titleForPath(pathname: string): string {
  return NAV_ITEMS.find((item) => item.href === pathname)?.label ?? "HelixScope";
}
