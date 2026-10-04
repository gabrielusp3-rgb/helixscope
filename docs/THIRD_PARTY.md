# Third-party components

Versions below are the ones measured in the tested HelixScope container.
License text for third-party software is not relicensed as HelixScope MIT.
Where this repository does not carry the upstream license file, the license
name is omitted rather than guessed.

| Component | Version | Role | Upstream |
| --- | --- | --- | --- |
| Python | 3.14.3 | Runtime | https://www.python.org/ |
| Streamlit | 1.58.0 | Interface | https://streamlit.io/ |
| Biopython | 1.87 | Sequence and protein calculations | https://biopython.org/ |
| Plotly | 6.9.0 | Charts | https://plotly.com/python/ |
| pandas | 2.3.3 | Tables | https://pandas.pydata.org/ |
| NumPy | 2.5.3 | Arrays | https://numpy.org/ |
| pytest | 9.1.1 | Tests inside the image | https://pytest.org/ |
| Hypothesis | 6.168.3 | Property tests | https://hypothesis.readthedocs.io/ |
| FastTree | 2.1.11 | Approximate maximum-likelihood trees | http://www.microbesonline.org/fasttree/ |
| IQ-TREE | 3.1.3 | Maximum-likelihood trees. The Dockerfile records the downloaded release as GPL-2.0-or-later. | https://github.com/iqtree/iqtree3 |
| BLAST+ | 2.17.0+ | Local similarity search against the tiny fixtures in the image | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ |
| US-align | 20260328 | Structure alignment. Citation recorded in `modules/usalign.py`. | https://zhanggroup.org/US-align/ |
| ViennaRNA | 2.7.2 | RNA folding through the Python module | https://www.tbi.univie.ac.at/RNA/ |
| NCBI E-utilities, ClinVar | remote | Record retrieval when the user asks | https://www.ncbi.nlm.nih.gov/ |
| Ensembl VEP | remote | Variant consequences when the user asks | https://www.ensembl.org/ |
| RCSB PDB | remote | Structure files and the RCSB alignment service | https://www.rcsb.org/ |
| EMBL-EBI | remote | Clustal Omega jobs when that MSA backend is selected | https://www.ebi.ac.uk/ |

DSSP, STRIDE, MSMS, EDTSurf and Cas-OFFinder are not installed in the tested
image. HelixScope reports them as not installed.
