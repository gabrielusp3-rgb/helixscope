# Optional scientific engines

HelixScope does not bundle IQ-TREE, mkdssp, FastTree, STRIDE, BLAST+ or ViennaRNA.
Place official binaries under this directory or set environment variables.
There is no recursive scan of `C:\`.

## IQ-TREE 3 (official)

Source: https://github.com/iqtree/iqtree3/releases  
Windows asset: `iqtree-3.1.3-Windows.zip`  
SHA-256: `6e0e11c7ed6197a0f4ebfbadafc7e82bebdf5b6e510c8a13f36183fe8eae3aef`  
License: GPL-2.0-or-later

Extract so the executable is one of:

- `tools/iqtree/bin/iqtree3.exe`
- `tools/iqtree-3.1.3-Windows/bin/iqtree3.exe`

Or set `HELIXSCOPE_IQTREE` to the full path of `iqtree3.exe`.

## mkdssp / DSSP 4 (official)

Upstream: https://github.com/PDB-REDO/dssp (BSD-2-Clause)
Local Windows build requires CMake and a C++20 compiler, which this workstation
does not provide. HelixScope therefore uses the official PDB-REDO DSSP HTTP API
(`https://pdb-redo.eu/dssp/do`, multipart POST) as a **remote** backend.

That path is `REMOTE_VALIDATED`, never `LIVE_VALIDATED` local mkdssp.

## EDTSurf (official Zhang lab)

Source: https://zhanggroup.org/EDTSurf/
Windows executable: `https://zhanggroup.org/EDTSurf/EDTSurf.exe`
Authors do not publish a checksum; HelixScope stores a local SHA-256 in
`tools/edtsurf/ORIGIN.json` after TLS retrieval and PE-header check.
License: permissive Zhang-lab disclaimer (retain copyright and citation).
Place `edtsurf.exe` at `tools/edtsurf/edtsurf.exe` or set `HELIXSCOPE_SURFACE`.
EDTSurf meshes are VDW/SAS/MS/SES triangulations, not Shrake-Rupley SASA.

## FastTree (official)

Source: https://morgannprice.github.io/fasttree/
Windows: `FastTree.exe` from that page (2.2.0, AVX2 required). Place at `tools/fasttree/FastTree.exe`
or set `HELIXSCOPE_FASTTREE`. Labeled approximate maximum likelihood, not IQ-TREE.

## US-align (official)

Source: https://zhanggroup.org/US-align/
Windows zip: https://zhanggroup.org/US-align/bin/module/USalignWin64.zip
Extract so `USalign.exe` is under `tools/usalign/` (depth <= 3) or set `HELIXSCOPE_USALIGN`.
The Win64 package observed here reports Version 20241108. Nucleic-acid comparison uses `-mol RNA`.

## BLAST+ (official NCBI)

Source: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
Windows tarball: `ncbi-blast-2.17.0+-x64-win64.tar.gz`

Local search runs when `HELIXSCOPE_BLAST_DB` points to a prefix with
`.nin`/`.nsq` (or protein equivalents), or when the declared fixtures
`tools/blast_db/helixscope_tiny_nucl` and `helixscope_tiny_prot` exist.
HelixScope never downloads nt/nr.

## ViennaRNA Python (official PyPI)

`pip install ViennaRNA==2.7.2` provides the `RNA` module (cp314 win_amd64 wheel).
That is not the `RNAfold` CLI. The CLI stays NOT_INSTALLED unless `RNAfold.exe` is present.

## STRIDE

Upstream: Heinig STRIDE (academic source). HelixScope has no official Windows
executable on this host. Absence is NOT_INSTALLED. Do not substitute DSSP.
Set `HELIXSCOPE_STRIDE` only if you have a legal STRIDE binary.

## MSMS (official Scripps CCSB, academic/teaching)

Source: https://ccsb.scripps.edu/msms/downloads/
Windows zip: `msms_win32_2.6.1.zip`
License: free for academic/teaching use; commercial research needs an agreement
with Dr Sanner. HelixScope does not commit `msms.exe`.

Extract so `msms.exe` is at `tools/msms/msms.exe` (depth <= 3) or set
`HELIXSCOPE_MSMS`. HelixScope writes xyzr with Bondi radii (not the MSMS
`pdb_to_xyzr` tables) and parses `.vert`/`.face`. SES from MSMS is not
Shrake-Rupley SASA and not EDTSurf. When both MSMS and EDTSurf are present,
EDTSurf remains the default surface engine.

## ViennaRNA CLI (official TBI Windows installer)

Source: https://www.tbi.univie.ac.at/RNA/
Windows installer: ViennaRNA Package v2.7.2 x86_64 (same version as the
installed Python wheel). The installer does not add PATH.

Place `RNAfold.exe` at `tools/viennarna/RNAfold.exe` together with the MinGW
DLLs it imports (`libstdc++-6.dll`, `libgcc_s_seh-1.dll`,
`libwinpthread-1.dll`), or under `C:\Program Files\ViennaRNA\RNAfold.exe`,
or set `HELIXSCOPE_RNAFOLD`.
The official installer may remain beside those files; HelixScope can unpack
the NSIS streams without running the vendor GUI.
ViennaRNA Python (`import RNA`) is a separate backend.

## Cas-OFFinder 2.4.1 (official)

Source: https://github.com/snugel/cas-offinder/releases/tag/2.4.1  
Windows asset: `cas-offinder_windows_x86-64.zip`  
License: BSD (upstream). CLI: `cas-offinder {input} {C|G|A} {output}`

HelixScope prefers device `C` (CPU) when a CPU OpenCL ICD is present. If only a
GPU OpenCL device exists, device `G` is used and labelled as GPU. GPU is not a
requirement. Cas-OFFinder 3 is not used. Native DNA/RNA bulges are UNAVAILABLE
on 2.4.1. HelixScope does not reimplement Cas-OFFinder.

Extract so the executable is one of:

- `tools/cas-offinder/cas-offinder.exe`
- `tools/cas-offinder.exe`

Or set `HELIXSCOPE_CAS_OFFINDER` to the full path of `cas-offinder.exe`.
The zip SHA-256 is recorded in `modules/genome_download.py` (`CAS_OFFINDER_ZIP_SHA256`)
after hashing the official GitHub asset; empty SHA means install is refused.

## Environment

Copy `.env.example` values into your shell. Do not commit `.env` or machine paths.
