# Data sources

HelixScope uses a remote source only for the request you submit. Local
sequence statistics do not contact these services.

## NCBI

- Role: nucleotide, protein, gene and PubMed records; ClinVar via E-utilities; remote BLAST.
- Access: HTTPS to NCBI when you search, fetch or submit BLAST. A contact email is required for Entrez.
- Optional key: `NCBI_API_KEY` or `ENTREZ_API_KEY`. Without a key, HelixScope stays inside the public request rate.
- Limit: an empty result or a network error stays empty or an error. It is not rewritten as a computed sequence.
- URL: https://www.ncbi.nlm.nih.gov/

The container does not include NCBI nt or nr.

## Ensembl

- Role: VEP consequence for a variant you choose to annotate.
- Access: HTTPS on demand.
- Limit: a VEP consequence is `RETRIEVED`. It is not a pathogenicity call.
- URL: https://www.ensembl.org/

## EMBL-EBI

- Role: Clustal Omega MSA (`ebi_clustalo`) and InterPro domain records.
- Access: HTTPS on demand.
- Limit: the tested container has no local MAFFT, MUSCLE or Clustal Omega binary. Offline MSA submission to EBI fails as a network error.
- URL: https://www.ebi.ac.uk/

## UniProt

- Role: protein cross-references used with domain annotation.
- Access: HTTPS on demand.
- URL: https://www.uniprot.org/

## RCSB PDB

- Role: experimental mmCIF files, entry metadata, and the RCSB Alignment API.
- Access: HTTPS on demand.
- Limit: an RCSB file is retrieved experimental data when the entry is experimental. It is not recalculated by HelixScope.
- URL: https://www.rcsb.org/

## PDB-REDO

- Role: remote DSSP (`mkdssp`) for secondary-structure assignment of coordinates you send.
- Access: HTTPS to `https://pdb-redo.eu/dssp/`.
- Limit: this is not a local DSSP binary. If the service fails, the assignment is unavailable.

## AlphaFold DB

- Role: predicted structure files requested by a UniProt accession.
- Access: HTTPS on demand.
- Limit: status is predicted, not experimental.
- URL: https://alphafold.ebi.ac.uk/

## Not a data source

FastTree, IQ-TREE, BLAST+ and US-align run locally on files already in the
session or in the image. They are engines, not databases.
