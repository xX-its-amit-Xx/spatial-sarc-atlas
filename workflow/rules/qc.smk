# Spot-level QC and filtering.

rule qc_visium:
    input:
        sentinel = "data/raw/{sample}/.fetched",
    output:
        adata    = "results/visium/{sample}/adata_qc.h5ad",
        report   = "results/figures/{sample}_qc.pdf",
    log:
        "logs/qc/{sample}.log",
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p results/visium/{wildcards.sample} results/figures logs/qc
        python -c "
import scanpy as sc, matplotlib.pyplot as plt
from pathlib import Path
import spatialdata_io as sdio

raw = Path('data/raw/{wildcards.sample}')
# Heuristic: try the standard 10x Space Ranger output layout first.
try:
    adata = sc.read_visium(raw)
except Exception:
    # Fall back to spatialdata-io for non-standard layouts.
    sdata = sdio.visium(raw)
    adata = sdata.tables['table']

adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                           log1p=False, inplace=True)

# Filter spots and genes with conservative defaults.
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs['pct_counts_mt'] < 25].copy()

fig, axes = plt.subplots(1, 3, figsize=(12, 3))
axes[0].hist(adata.obs['total_counts'], bins=60); axes[0].set_title('UMI / spot')
axes[1].hist(adata.obs['n_genes_by_counts'], bins=60); axes[1].set_title('Genes / spot')
axes[2].hist(adata.obs['pct_counts_mt'], bins=60); axes[2].set_title('%MT')
fig.tight_layout()
fig.savefig('{output.report}', bbox_inches='tight')

adata.write('{output.adata}')
" > {log} 2>&1
        """
