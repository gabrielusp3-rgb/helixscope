# Genome references (local)

HelixScope does not ship GRCh38, T2T-CHM13 or GRCm39. This directory holds
checksum-verified FASTA files after an explicit Download action in the CRISPR
tab (or `genome_download.download_and_install_assembly`).

READY means:

1. the compressed `*_genomic.fna.gz` MD5 matched NCBI `md5checksums.txt`;
2. gzip was decompressed under a 20 GiB cap;
3. a samtools-style FAI index and `.fai.meta.json` (FASTA SHA-256) were written.

A `.partial` file is never READY. Opening the CRISPR tab does not download
anything. Do not commit FASTA, `.gz`, `.fai` or `manifest.json` from this tree.

Set `HELIXSCOPE_REFERENCE_DIR` to use another location. Disk: budget several
gigabytes per assembly (compressed ~0.8-1 GB plus uncompressed FASTA).

NCBI RefSeq data: cite the assembly accession (GCF_000001405.40, GCF_009914755.1,
GCF_000001635.27). A file named `hg38.fa` is not proof of GRCh38.
