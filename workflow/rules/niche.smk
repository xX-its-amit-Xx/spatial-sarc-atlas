# Niche identification and tumor-stroma interface analysis.

rule niche_analysis:
    input:
        adata    = "results/visium/{sample}/adata_c2l.h5ad",
    output:
        adata    = "results/visium/{sample}/niche_metrics.h5ad",
        nhood    = "results/figures/{sample}_nhood_enrichment.pdf",
        cooc     = "results/figures/{sample}_cooccurrence.pdf",
    log:
        "logs/niche/{sample}.log",
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p logs/niche
        python -c "
import scanpy as sc, squidpy as sq, numpy as np
import matplotlib.pyplot as plt

adata = sc.read_h5ad('{input.adata}')
# Use the cell2location dominant-celltype call as the cluster_key.
if 'dominant_celltype' not in adata.obs:
    prop = adata.obsm.get('q05_cell_abundance_w_sf')
    if prop is None:
        raise SystemExit('cell2location output missing; run deconvolution first.')
    adata.obs['dominant_celltype'] = prop.idxmax(axis=1).astype('category')

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
sq.gr.nhood_enrichment(adata, cluster_key='dominant_celltype')
sq.gr.co_occurrence(adata, cluster_key='dominant_celltype')

fig1 = plt.figure(figsize=(6, 6))
sq.pl.nhood_enrichment(adata, cluster_key='dominant_celltype', method='ward',
                       ax=plt.gca(), title='Neighborhood enrichment')
fig1.savefig('{output.nhood}', bbox_inches='tight')

fig2 = plt.figure(figsize=(8, 4))
sq.pl.co_occurrence(adata, cluster_key='dominant_celltype', ax=plt.gca())
fig2.savefig('{output.cooc}', bbox_inches='tight')

adata.write('{output.adata}')
" > {log} 2>&1
        """
