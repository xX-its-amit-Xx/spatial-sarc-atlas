# spatial-sarc-atlas

A reproducible reanalysis pipeline for publicly available 10x Visium and NanoString CosMx spatial transcriptomics data on **sarcoma** samples.

## Why this exists

Most spatial-omics tutorials live in brain, breast, or colorectal cancer. The sarcoma tumor microenvironment (TME) — a histologically and molecularly heterogeneous group of mesenchymal malignancies — is comparatively understudied with spatial methods, even though stromal architecture is arguably more interpretive in sarcoma than in any other tumor class. This repo reanalyzes a handful of public sarcoma Visium and CosMx datasets end-to-end (QC → normalization → deconvolution → niche identification → tumor–stroma interface differential expression) so the methods are easy to reuse and the figures are easy to scroll. It is intentionally notebook-first.

## Quick start

### Conda
```bash
git clone https://github.com/xX-its-amit-Xx/spatial-sarc-atlas.git
cd spatial-sarc-atlas
conda env create -f environment.yml
conda activate spatial-sarc-atlas
snakemake --use-conda --cores 4 -p
```

### Docker
```bash
docker build -t spatial-sarc-atlas .
docker run --rm -it -v "$PWD":/work -w /work spatial-sarc-atlas \
    snakemake --cores 4 -p
```

### Just the notebooks
```bash
jupyter lab notebooks/
```

## Notebook gallery

Rendered with `jupyter-book` and deployed to GitHub Pages on every push to `main`:
**https://xX-its-amit-Xx.github.io/spatial-sarc-atlas/**

| # | Notebook | What it covers |
|---|---|---|
| 01 | [QC & filtering](notebooks/01_qc_and_filtering.ipynb) | Spot-level QC, mitochondrial content, tissue masking |
| 02 | [Normalization & clustering](notebooks/02_normalization_and_clustering.ipynb) | SCT vs. log-norm, Leiden, spatially aware clustering |
| 03 | [Cell-type deconvolution](notebooks/03_celltype_deconvolution.ipynb) | cell2location **and** Tangram, side-by-side |
| 04 | [Tumor–stroma interface](notebooks/04_tumor_stroma_interface.ipynb) | Boundary spots, DE at the interface |
| 05 | [Niche identification](notebooks/05_niche_identification.ipynb) | Squidpy `nhood_enrichment`, `co_occurrence`, niches |

## Datasets reanalyzed

Accessions are pinned in [config/samples.yaml](config/samples.yaml). All datasets are publicly accessible; this repo redistributes **no** raw data.

| Modality | Tumor type | Accession | Verification | Reference |
|---|---|---|---|---|
| 10x Visium | Angiosarcoma | GSE227469 | ✅ Verified 2026-05-20 (GPL24676) | Li et al., *Commun Biol* (2023) |
| CosMx WTX | Soft-tissue sarcoma demo | NanoString public release | Public download URL (no platform check) | NanoString CosMx data portal |
| scRNA reference | Synovial sarcoma | GSE131309 | Manual conversion to h5ad required | Jerby-Arnon et al., *Nat Med* (2021) |
| scRNA reference | Sarcoma (archival, multi-subtype) | dbGaP / GEO (see config) | Manual conversion to h5ad required | Subramanian et al., *Clin Cancer Res* (2024) |

> **Verification.** `scripts/fetch_public_data.py --verify-only` pings GEO and confirms each Visium accession advertises a Visium platform tag (GPL24676 / GPL30178 / GPL30172). The two original placeholder accessions (GSE190294, GSE205276) were rejected on 2026-05-20 — see the commented-out block in [config/samples.yaml](config/samples.yaml) for the rejection details and how to add new verified accessions.

## Repository layout

```
spatial-sarc-atlas/
├── README.md
├── LICENSE                  # MIT
├── environment.yml
├── Dockerfile
├── Snakefile
├── config/
│   ├── samples.yaml         # Public dataset accessions (GSE, Zenodo DOIs)
│   └── references.yaml      # Cell-type reference dataset locations
├── workflow/
│   └── rules/
│       ├── download.smk
│       ├── qc.smk
│       ├── normalize.smk
│       ├── deconvolution.smk
│       └── niche.smk
├── notebooks/               # The public artifact
├── scripts/
│   ├── fetch_public_data.py
│   ├── build_reference.py
│   ├── run_cell2location.py
│   └── run_tangram.py
├── data/
│   └── references/          # Not committed; built/downloaded on demand
├── results/
│   ├── figures/
│   └── reports/             # Per-sample auto-generated markdown summaries
└── docs/                    # jupyter-book source for GitHub Pages
```

## Cookbook — real-world recipes

Each recipe is a self-contained snippet that uses one of the reanalyzed datasets. Treat these as "common-task shortcuts" once you have run the pipeline once.

### Recipe 1 — Load a processed Visium AnnData and re-plot a marker gene

```python
import scanpy as sc

adata = sc.read_h5ad("results/visium/GSE227469_sample01/adata_normalized.h5ad")
sc.pl.spatial(adata, color=["VWF", "PECAM1", "MKI67"], size=1.4,
              ncols=3, cmap="magma", save="_GSE227469_endothelial.pdf")
```

*Use when:* you want to inspect endothelial markers in the angiosarcoma sample without rerunning QC.

### Recipe 2 — Get cell2location proportions as a tidy data frame

```python
import scanpy as sc, pandas as pd

adata = sc.read_h5ad("results/visium/GSE227469_sample01/adata_c2l.h5ad")
prop = adata.obsm["q05_cell_abundance_w_sf"].copy()
prop.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in prop.columns]
prop = prop.div(prop.sum(axis=1), axis=0)
prop.index.name = "spot"
prop.reset_index().to_csv("results/reports/GSE227469_sample01_proportions.csv", index=False)
```

*Use when:* you want to hand a downstream collaborator a CSV of estimated cell-type fractions per spot.

### Recipe 3 — Identify tumor–stroma interface spots in one pass

```python
import scanpy as sc, squidpy as sq, numpy as np

adata = sc.read_h5ad("results/visium/GSE227469_sample01/adata_c2l.h5ad")
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

tumor_frac  = adata.obs["fraction_tumor"]
stromal_frac = adata.obs["fraction_caf"] + adata.obs["fraction_endothelial"]

interface = (tumor_frac.between(0.25, 0.75)) & (stromal_frac.between(0.25, 0.75))
adata.obs["compartment"] = np.where(
    interface, "interface",
    np.where(tumor_frac > 0.75, "tumor",
    np.where(stromal_frac > 0.75, "stroma", "other"))
)

sc.tl.rank_genes_groups(adata, "compartment", method="wilcoxon",
                        groups=["interface"], reference="tumor")
sc.pl.rank_genes_groups(adata, n_genes=20, save="_interface_vs_tumor.pdf")
```

*Use when:* you want a quick boundary-DE result without opening Notebook 04.

### Recipe 4 — Squidpy neighborhood enrichment on existing clusters

```python
import scanpy as sc, squidpy as sq

adata = sc.read_h5ad("results/visium/GSE227469_sample01/adata_normalized.h5ad")
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(adata, cluster_key="leiden",
                       method="ward", figsize=(6, 6),
                       save="GSE227469_sample01_nhood.pdf")
```

*Use when:* you want the niche-enrichment matrix without rerunning the full Notebook 05.

### Recipe 5 — Compare cell2location vs. Tangram for one cell type

```python
import scanpy as sc, matplotlib.pyplot as plt

c2l = sc.read_h5ad("results/visium/GSE227469_sample01/adata_c2l.h5ad")
tg  = sc.read_h5ad("results/visium/GSE227469_sample01/adata_tangram.h5ad")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sc.pl.spatial(c2l, color="fraction_caf", ax=axes[0], show=False,
              title="cell2location — CAF", cmap="viridis")
sc.pl.spatial(tg,  color="CAF", ax=axes[1], show=False,
              title="Tangram — CAF", cmap="viridis")
fig.savefig("results/figures/GSE227469_sample01_caf_method_comparison.pdf",
            bbox_inches="tight")
```

*Use when:* you want a quick visual sanity-check that the two deconvolvers agree where it matters.

## Reproducibility notes

- `environment.yml` pins major versions; `Dockerfile` pins the conda env at build time.
- Each Snakemake rule re-uses the conda env, so the pipeline works on HPC with `--profile slurm`.
- GPU is **optional** but recommended for cell2location and scVI training. The pipeline falls back to CPU.
- Per-sample auto-generated `results/reports/<sample>_summary.md` records spot counts, median UMI, top niche, and the package versions used.

## Limitations

This repo is honest about what it does *not* do yet:

- **No MERFISH or Xenium yet.** Only Visium (FFPE + FF) and CosMx are supported. Adding `spatialdata-io` parsers for those is a planned next step.
- **Reference atlases are sarcoma-leaning, not sarcoma-perfect.** Synovial sarcoma and angiosarcoma are well-covered; rare subtypes (e.g., epithelioid sarcoma, alveolar soft-part sarcoma) inherit cell-state labels from the closest available subtype.
- **No batch correction across studies by default.** The notebooks treat each sample independently. A multi-sample integration step (Harmony/scVI) is shown but not part of the auto pipeline.
- **Niche labels are unsupervised.** They are stable but not biologically annotated automatically — interpretation is a manual step in Notebook 05.
- **No clinical metadata harmonization.** Survival and treatment annotations from the original studies are not normalized into a common schema.

## Acknowledgments and citation

If you use the figures, code, or interpretations here, please cite the original studies whose data is reanalyzed (see [config/samples.yaml](config/samples.yaml) for citations alongside each accession) **and** the methods packages: scanpy, squidpy, cell2location, Tangram, scvi-tools.

This repo itself can be cited as:

> Shenoy A. *spatial-sarc-atlas: a reproducible reanalysis pipeline for sarcoma spatial transcriptomics.* 2026. https://github.com/xX-its-amit-Xx/spatial-sarc-atlas

## License

MIT — see [LICENSE](LICENSE).
