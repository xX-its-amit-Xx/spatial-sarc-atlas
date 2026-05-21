# Normalization, HVG selection, dim reduction, clustering.

rule normalize_visium:
    input:
        adata    = "results/visium/{sample}/adata_qc.h5ad",
    output:
        adata    = "results/visium/{sample}/adata_normalized.h5ad",
    log:
        "logs/normalize/{sample}.log",
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p logs/normalize
        python -c "
import scanpy as sc

adata = sc.read_h5ad('{input.adata}')
adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3',
                            layer='counts')
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.6, key_added='leiden')

adata.write('{output.adata}')
" > {log} 2>&1
        """
