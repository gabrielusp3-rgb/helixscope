# Limitations

## Local

- The technical DNA/RNA/protein input cap is 200,000 residues after cleaning,
  and 250,000 raw characters before that.
- Pairwise alignment stops at 8,000 residues per sequence. Dot plots stop at
  2,000.
- Global RNA fold stops at 600 nt.
- The illustrative 3D helix stops at 400 nt and is not the analytical result.
- History is kept in the browser session, up to 100 analyses. It is not a
  shared database.
- SEG low-complexity masking is not implemented. Shannon entropy is a different
  number.

## Docker

- The tested Compose limits are 64 processes and 1 GiB.
- The root filesystem is read-only only when Compose (or an equivalent run)
  sets that flag. The image alone does not enforce it.
- Cas-OFFinder, local DSSP, STRIDE, MSMS and EDTSurf are absent.
- BLAST databases inside the image are the tiny fixtures only.

## Remote APIs

- NCBI, Ensembl, EMBL-EBI, UniProt, RCSB, PDB-REDO and AlphaFold DB require
  outbound HTTPS and can change or refuse a request.
- Without `NCBI_API_KEY`, Entrez uses the public rate limit.
- An MSA cannot be computed offline in the tested image, because the only
  offered backend is EBI Clustal Omega.

## Engines

- FastTree and IQ-TREE produce a tree for the alignment and model you gave
  them. They do not produce bootstrap values unless that option was run.
- IQ-TREE does not accept p-distance.
- Local `blastn` uses DUST. Low-complexity queries can be `NO_HITS`.
- US-align in this interface is monomer alignment (`-mm 0`).
- OpenCL and a GPU are not available in the container.

## Platforms

- Streamlit Community Cloud does not build this Dockerfile. It also cannot
  pin Python to 3.14.3 from the repository alone.
- A Render free instance is 512 MB and 0.1 CPU, and it sleeps after 15 minutes
  without traffic. That is not the 1 GiB container used for the tests.
- No public URL is provided by this repository.

## Scientific

- No pathogenicity score.
- No genome-wide CRISPR specificity.
- No claim that a tree is the true history of the sequences.
- No claim that an RNA secondary structure is an experimental 3D structure.
- No claim that a 1CRN self-alignment describes any other pair.
