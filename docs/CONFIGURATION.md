# Configuration

Only variables the code reads are listed. None of them has a default secret.

| Variable | Used for | If unset |
| --- | --- | --- |
| `PORT` | TCP port inside the container | 8501 |
| `NCBI_API_KEY` | Optional NCBI E-utilities key | Public rate limit |
| `ENTREZ_API_KEY` | Accepted as another name for the NCBI key | Public rate limit |
| `HELIXSCOPE_REFERENCE_DIR` | Local reference-genome directory | `data/references` under the checkout |
| `HELIXSCOPE_JOBS_DIR` | Local genome-job directory | `data/jobs` |
| `HELIXSCOPE_STRUCTURE_CACHE` | Cache directory for retrieved structure files | The module default under the checkout |
| `HELIXSCOPE_BLAST_DB` | Prefix of a local BLAST database | The tiny fixture prefix inside the image, when that file exists |
| `HELIXSCOPE_TOOLS_DIR` | Optional directory of local executables | `PATH` only |
| `HELIXSCOPE_FASTTREE` | Explicit FastTree path | `PATH` |
| `HELIXSCOPE_IQTREE` | Explicit IQ-TREE path | `PATH` |
| `HELIXSCOPE_RNAFOLD` | Explicit RNAfold executable, distinct from the ViennaRNA Python module | Not required for the Python folder |

`.env.example` lists the optional engine paths without values. Do not commit
a filled `.env`.

`PORT` is read by `docker/entrypoint.py` and `docker/healthcheck.py`. The
Streamlit process on your own machine, started without Docker, uses
`--server.port 8501` unless you pass another flag. It does not read `PORT`.
