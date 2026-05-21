# spatial-sarc-atlas — top-level Snakefile
#
# Run with:
#   snakemake --use-conda --cores 4 -p
#
# Important targets:
#   snakemake all              -> everything (default)
#   snakemake fetch            -> download raw data only
#   snakemake reference        -> build sarcoma scRNA reference only
#   snakemake reports          -> render per-sample summary markdowns
#   snakemake notebooks        -> execute all notebooks via papermill
#
# Pipeline is notebook-aware: every per-sample step writes an h5ad that the
# notebooks pick up. This means you can run the pipeline end-to-end OR open
# the notebooks and run interactively from any checkpoint.

from pathlib import Path
import yaml

configfile: "config/samples.yaml"

# ---------------------------------------------------------------------------
# Sample expansion
# ---------------------------------------------------------------------------
VISIUM = [s["id"] for s in config.get("visium", [])]
COSMX  = [s["id"] for s in config.get("cosmx", [])]
ALL_SAMPLES = VISIUM + COSMX

RESULTS = Path("results")
DATA    = Path("data")

# ---------------------------------------------------------------------------
# Top-level targets
# ---------------------------------------------------------------------------
rule all:
    input:
        # Per-sample deconvolved AnnData (both methods)
        expand(str(RESULTS / "visium/{sample}/adata_c2l.h5ad"),     sample=VISIUM),
        expand(str(RESULTS / "visium/{sample}/adata_tangram.h5ad"), sample=VISIUM),
        # Per-sample niche analyses
        expand(str(RESULTS / "visium/{sample}/niche_metrics.h5ad"), sample=VISIUM),
        # Per-sample report
        expand(str(RESULTS / "reports/{sample}_summary.md"),        sample=VISIUM),

rule fetch:
    input:
        expand(str(DATA / "raw/{sample}/.fetched"), sample=ALL_SAMPLES),

rule reference:
    input:
        "data/references/sarcoma_reference_signatures.h5ad",

rule notebooks:
    input:
        expand("results/notebooks/{nb}.ipynb",
               nb=["01_qc_and_filtering",
                   "02_normalization_and_clustering",
                   "03_celltype_deconvolution",
                   "04_tumor_stroma_interface",
                   "05_niche_identification"])

rule reports:
    input:
        expand(str(RESULTS / "reports/{sample}_summary.md"), sample=VISIUM),

# ---------------------------------------------------------------------------
# Per-stage rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/download.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/normalize.smk"
include: "workflow/rules/deconvolution.smk"
include: "workflow/rules/niche.smk"

# ---------------------------------------------------------------------------
# Notebook execution (papermill)
# ---------------------------------------------------------------------------
rule run_notebook:
    input:
        nb="notebooks/{nb}.ipynb",
    output:
        nb="results/notebooks/{nb}.ipynb",
    log:
        "logs/notebooks/{nb}.log",
    conda:
        "environment.yml"
    shell:
        """
        mkdir -p results/notebooks logs/notebooks
        papermill {input.nb} {output.nb} \
            -p results_dir results \
            -p data_dir data \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------------
# Per-sample auto-summary
# ---------------------------------------------------------------------------
rule sample_summary:
    input:
        adata   = RESULTS / "visium/{sample}/adata_c2l.h5ad",
        niche   = RESULTS / "visium/{sample}/niche_metrics.h5ad",
    output:
        report  = RESULTS / "reports/{sample}_summary.md",
    log:
        "logs/reports/{sample}.log",
    conda:
        "environment.yml"
    script:
        "scripts/write_sample_summary.py"
