# Download raw public data for each sample.

rule fetch_sample:
    output:
        sentinel = touch("data/raw/{sample}/.fetched"),
    params:
        config   = "config/samples.yaml",
    log:
        "logs/fetch/{sample}.log",
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p data/raw/{wildcards.sample} logs/fetch
        python scripts/fetch_public_data.py \
            --config {params.config} \
            --sample {wildcards.sample} \
            --out data/raw/{wildcards.sample} \
            > {log} 2>&1
        """

rule build_reference:
    output:
        combined   = "data/references/sarcoma_reference_combined.h5ad",
        signatures = "data/references/sarcoma_reference_signatures.h5ad",
    params:
        config = "config/references.yaml",
    log:
        "logs/reference/build.log",
    conda:
        "../../environment.yml"
    shell:
        """
        mkdir -p data/references logs/reference
        python scripts/build_reference.py \
            --config {params.config} \
            --out-combined {output.combined} \
            --out-signatures {output.signatures} \
            > {log} 2>&1
        """
