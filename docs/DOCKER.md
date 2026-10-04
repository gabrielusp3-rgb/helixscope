# Docker

The maintained container is a local Linux image. Building it does not publish
it and does not create a hosting account.

## Tested image

- Base tag: `python:3.14.3-slim-trixie`
- Base digest recorded in the Dockerfile:
  `sha256:5e59aae31ff0e87511226be8e2b94d78c58f05216efda3b07dbbed938ec8583b`
- Python 3.14.3, Streamlit 1.58.0, Biopython 1.87
- Plotly 6.9.0, pandas 2.3.3, NumPy 2.5.3
- FastTree 2.1.11, IQ-TREE 3.1.3, BLAST+ 2.17.0+, US-align 20260328,
  ViennaRNA 2.7.2
- User `helix`, uid 10001
- Platform `linux/amd64`

The image id changes when the Dockerfile changes. The digest above is the
base image, not the HelixScope image id.

## Run locally

```text
docker compose up --build
```

Open http://127.0.0.1:8501 or http://localhost:8501.

Compose sets `PORT=8501` and publishes that port only on the host loopback.
The process binds `0.0.0.0` inside the container.

Health, for that Compose port:

```text
http://127.0.0.1:8501/_stcore/health
```

The body is `ok` when Streamlit is up.

## PORT

`docker/entrypoint.py` reads `PORT`.

- Unset: 8501
- `8501`, `9000`, `10000`: that port
- Empty, non-numeric, or outside 1..65535: startup exits with
  `PORT must be an integer from 1 to 65535.`

The health check calls `http://127.0.0.1:$PORT/_stcore/health`. A second
container can use another port without rebuilding, for example:

```text
docker run --rm -e PORT=9000 -p 127.0.0.1:9000:9000 helixscope:local
```

Compose itself stays on 8501. Changing only `PORT` in Compose without changing
the published port mapping will make the host URL miss the process.

## Runtime limits in Compose

- Read-only root filesystem
- Writable temporary mounts for `/tmp` and `/home/helix`
- 64 processes
- 1 GiB memory
- No extra capabilities and no privileged mode

A 512 MB free hosting plan is smaller than this tested limit. It was not
validated as a substitute.

## What the image does not contain

Cas-OFFinder, local DSSP, STRIDE, MSMS, EDTSurf, NCBI nt, NCBI nr, and
reference assemblies. Secrets are not baked into the image.
