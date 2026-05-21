# Cell-type deconvolution with two methods so users can compare.

rule cell2location:
    input:
        adata      = "results/visium/{sample}/adata_normalized.h5ad",
        signatures = "data/references/sarcoma_reference_signatures.h5ad",
    output:
        adata      = "results/visium/{sample}/adata_c2l.h5ad",
        figure     = "results/figures/{sample}_c2l_proportions.pdf",
    log:
        "logs/c2l/{sample}.log",
    threads: 4
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p logs/c2l
        python scripts/run_cell2location.py \
            --input {input.adata} \
            --signatures {input.signatures} \
            --output {output.adata} \
            --figure {output.figure} \
            > {log} 2>&1
        """

rule tangram:
    input:
        adata      = "results/visium/{sample}/adata_normalized.h5ad",
        reference  = "data/references/sarcoma_reference_combined.h5ad",
    output:
        adata      = "results/visium/{sample}/adata_tangram.h5ad",
        figure     = "results/figures/{sample}_tangram_proportions.pdf",
    log:
        "logs/tangram/{sample}.log",
    threads: 4
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p logs/tangram
        python scripts/run_tangram.py \
            --input {input.adata} \
            --reference {input.reference} \
            --output {output.adata} \
            --figure {output.figure} \
            > {log} 2>&1
        """
