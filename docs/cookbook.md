# Cookbook

Real-world recipes that use the outputs of the main pipeline. Each recipe assumes you have run `snakemake --use-conda --cores 4` once and have processed AnnData files under `results/visium/<sample>/`.

## Recipe 1 — Replot a marker gene without rerunning QC

```python
import scanpy as sc

adata = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_normalized.h5ad")
sc.pl.spatial(adata, color=["VWF", "PECAM1", "MKI67"], size=1.4,
              ncols=3, cmap="magma", save="_GSE227469_endothelial.pdf")
```

## Recipe 2 — Export cell-type fractions as CSV

```python
import scanpy as sc

adata = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_c2l.h5ad")
prop = adata.obsm["q05_cell_abundance_w_sf"].copy()
prop.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in prop.columns]
prop = prop.div(prop.sum(axis=1), axis=0)
prop.to_csv("results/reports/GSE227469_angiosarcoma_01_proportions.csv")
```

## Recipe 3 — Identify tumor–stroma interface spots inline

```python
import scanpy as sc
import numpy as np

adata = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_c2l.h5ad")
tumor  = adata.obs["fraction_malignant"]
stroma = adata.obs["fraction_caf"] + adata.obs["fraction_endothelial"]
adata.obs["compartment"] = np.where(
    (tumor.between(0.25, 0.75)) & (stroma.between(0.25, 0.75)), "interface",
    np.where(tumor > 0.75, "tumor",
    np.where(stroma > 0.75, "stroma", "other"))
)
sc.tl.rank_genes_groups(adata, "compartment", method="wilcoxon",
                        groups=["interface"], reference="tumor")
```

## Recipe 4 — Squidpy neighborhood enrichment

```python
import scanpy as sc, squidpy as sq

adata = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_normalized.h5ad")
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(adata, cluster_key="leiden", method="ward")
```

## Recipe 5 — Compare cell2location and Tangram for one cell type

```python
import scanpy as sc, matplotlib.pyplot as plt

c2l = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_c2l.h5ad")
tg  = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_tangram.h5ad")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sc.pl.spatial(c2l, color="fraction_caf", ax=axes[0], show=False,
              title="cell2location — CAF", cmap="viridis")
sc.pl.spatial(tg,  color="fraction_caf", ax=axes[1], show=False,
              title="Tangram — CAF", cmap="viridis")
```

## Recipe 6 — Boundary-aware DE on a custom compartment

```python
import scanpy as sc
adata = sc.read_h5ad("results/visium/GSE227469_angiosarcoma_01/adata_c2l.h5ad")
adata.obs["my_compartment"] = (adata.obs["fraction_macrophage"] > 0.2).map(
    {True: "tam_rich", False: "tam_low"}
)
sc.tl.rank_genes_groups(adata, "my_compartment", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=20)
```
