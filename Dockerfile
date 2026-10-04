# HelixScope runtime. Official Linux sources only. No Windows executables.
# FastTree source: http://www.microbesonline.org/fasttree/FastTree.c (FT_VERSION 2.1.11)
# IQ-TREE: https://github.com/iqtree/iqtree3/releases/tag/v3.1.3 GPL-2.0-or-later
# BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/
# US-align: https://zhanggroup.org/US-align/bin/module/USalignLinux64.zip
# ViennaRNA: PyPI ViennaRNA==2.7.2 cp314 manylinux wheel
# Base decision, 2026-10-03:
# python:3.14.3-slim-bookworm@sha256:c6f0b5b3a167963de3cc7cd97fe1a5d07105c8d27f47507ab31d648b533b05c7
# plus apt-get upgrade still left CRITICAL CVE-2026-13221 and CVE-2026-12087
# in perl-base 5.36.0-7+deb12u3, with no Bookworm fix.
# python:3.14.3-slim-trixie has those fixes in perl 5.40.1-6+deb13u1, and
# openssl/glibc security updates that clear the other criticals. Python stays 3.14.3.
# Tag python:3.14.3-slim-trixie, digest below.

FROM python:3.14.3-slim-trixie@sha256:5e59aae31ff0e87511226be8e2b94d78c58f05216efda3b07dbbed938ec8583b AS engines

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        unzip \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL -o /tmp/FastTree.c http://www.microbesonline.org/fasttree/FastTree.c \
    && echo "9026ae550307374be92913d3098f8d44187d30bea07902b9dcbfb123eaa2050f  /tmp/FastTree.c" | sha256sum -c - \
    && gcc -O3 -finline-functions -funroll-loops -o /opt/FastTree /tmp/FastTree.c -lm \
    && rm /tmp/FastTree.c

RUN curl -fsSL -o /tmp/iqtree.tar.gz \
        https://github.com/iqtree/iqtree3/releases/download/v3.1.3/iqtree-3.1.3-Linux-intel.tar.gz \
    && echo "ac87dee78d06b67a1be87fff4a325358d038b5ae947308e52b3cf23829521aa8  /tmp/iqtree.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/iqtree.tar.gz -C /tmp \
    && find /tmp -type f -name iqtree3 -executable -exec cp {} /opt/iqtree3 \; \
    && rm -rf /tmp/iqtree.tar.gz /tmp/iqtree-3.1.3-Linux-intel

RUN curl -fsSL -o /tmp/blast.tar.gz \
        https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz \
    && echo "bdec166721de3b55f90a3badc83538e8  /tmp/blast.tar.gz" | md5sum -c - \
    && tar -xzf /tmp/blast.tar.gz -C /tmp \
    && mkdir -p /opt/blast \
    && cp /tmp/ncbi-blast-2.17.0+/bin/blastn \
        /tmp/ncbi-blast-2.17.0+/bin/blastp \
        /tmp/ncbi-blast-2.17.0+/bin/blastx \
        /tmp/ncbi-blast-2.17.0+/bin/tblastn \
        /tmp/ncbi-blast-2.17.0+/bin/tblastx \
        /tmp/ncbi-blast-2.17.0+/bin/makeblastdb \
        /opt/blast/ \
    && rm -rf /tmp/blast.tar.gz /tmp/ncbi-blast-2.17.0+

RUN curl -fsSL -A "Mozilla/5.0 (compatible; HelixScope-build/1.0)" \
        -o /tmp/usalign.zip \
        https://zhanggroup.org/US-align/bin/module/USalignLinux64.zip \
    && echo "af5ad073c0f732de22e6dd62c5abf25b2cbad9c791d436fff32231c61576b7a3  /tmp/usalign.zip" | sha256sum -c - \
    && unzip -j /tmp/usalign.zip USalign/USalign -d /opt \
    && chmod 755 /opt/USalign \
    && rm /tmp/usalign.zip

FROM python:3.14.3-slim-trixie@sha256:5e59aae31ff0e87511226be8e2b94d78c58f05216efda3b07dbbed938ec8583b

RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && apt-get upgrade -y --no-install-recommends \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 10001 helix \
    && useradd --uid 10001 --gid 10001 --create-home --shell /usr/sbin/nologin helix

COPY --from=engines /opt/FastTree /opt/iqtree3 /opt/USalign /usr/local/bin/
COPY --from=engines /opt/blast/ /usr/local/bin/

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt \
    && pip install --no-cache-dir "ViennaRNA==2.7.2" \
    && python -m pip uninstall -y pip \
    && rm -rf /usr/local/lib/python3.14/site-packages/pip \
        /usr/local/lib/python3.14/site-packages/pip-*.dist-info

COPY app.py /app/app.py
COPY modules /app/modules
COPY ui /app/ui
COPY helixscope_core /app/helixscope_core
COPY .streamlit /app/.streamlit
COPY tests/fixtures/1CRN.cif tests/fixtures/1BNA.cif tests/fixtures/1RNA.cif /app/fixtures/
COPY docker /app/docker

RUN mkdir -p /app/tools/blast_db /app/data/references /app/data/jobs /home/helix \
    && printf '%s\n' '>tiny_ref' 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT' > /tmp/tiny_nucl.fa \
    && printf '%s\n' '>tiny_prot' 'ACDEFGHIKLMNPQRSTVWYACDEFGHIKL' > /tmp/tiny_prot.fa \
    && makeblastdb -in /tmp/tiny_nucl.fa -dbtype nucl -out /app/tools/blast_db/helixscope_tiny_nucl \
    && makeblastdb -in /tmp/tiny_prot.fa -dbtype prot -out /app/tools/blast_db/helixscope_tiny_prot \
    && rm -f /tmp/tiny_nucl.fa /tmp/tiny_prot.fa \
    && chown -R helix:helix /app/data /app/tools /home/helix \
    && chmod 755 /usr/local/bin/FastTree /usr/local/bin/iqtree3 /usr/local/bin/USalign \
        /usr/local/bin/blastn /usr/local/bin/blastp /usr/local/bin/blastx \
        /usr/local/bin/tblastn /usr/local/bin/tblastx /usr/local/bin/makeblastdb

USER helix

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    HOME=/home/helix \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_LOGGER_HIDE_WELCOME_MESSAGE=true

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=5s --start-period=40s --retries=3 \
    CMD ["python", "/app/docker/healthcheck.py"]

CMD ["python", "/app/docker/entrypoint.py"]
