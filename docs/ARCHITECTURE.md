# Architecture

```text
Browser
  -> Streamlit process (app.py, ui/)
       -> st.session_state     history and the active module, per browser session
       -> st.cache_data        pure calculations, keyed by their arguments
       -> modules/             science, parsers, provenance
            -> Python libraries (Biopython, ViennaRNA, Plotly)
            -> subprocess list  FastTree, IQ-TREE, BLAST+, US-align
            -> urllib allowlist NCBI, RCSB, PDB-REDO, EMBL-EBI, Ensembl, UniProt
```

Biology is not implemented in the browser. The frozen Next.js tree under
`apps/web` is not on this path.

## Processes

The container entrypoint replaces itself with Streamlit (`os.execv`). Engine
calls use an argument list. They do not pass the user sequence through a
shell string.

Temporary files for those engines live under the system temporary directory
and are removed by the calling module. Compose puts `/tmp` on a temporary
mount.

## Detection

On startup and on the engines screen, HelixScope asks `PATH` for known
basenames. A path containing `..` is rejected. A missing tool returns
`not_installed`.

## Remote calls

URL checks require HTTPS and a host plus path prefix that the module lists.
Link-local addresses, loopback, private ranges as literal URLs, and
`file:`, `ftp:` and `data:` URLs are not accepted by those allowlists.

## What is ephemeral

Session history, pasted sequences and unsaved uploads disappear when the
session ends. Reference genomes and job files, if you install any, use
`HELIXSCOPE_REFERENCE_DIR` and `HELIXSCOPE_JOBS_DIR`, or `data/references`
and `data/jobs`. Those directories are not part of the Git tree.
