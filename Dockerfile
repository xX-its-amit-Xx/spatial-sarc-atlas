FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="spatial-sarc-atlas"
LABEL org.opencontainers.image.description="Reproducible reanalysis pipeline for sarcoma spatial transcriptomics"
LABEL org.opencontainers.image.source="https://github.com/xX-its-amit-Xx/spatial-sarc-atlas"
LABEL org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        wget \
        ca-certificates \
        build-essential \
        procps \
        unzip \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
WORKDIR /work

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . /work/

CMD ["snakemake", "--cores", "4", "-p"]
