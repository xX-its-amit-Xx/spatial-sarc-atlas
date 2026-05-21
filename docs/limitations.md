# Limitations

This pipeline is intentionally scoped narrow. Things it does *not* do yet:

- **MERFISH / Xenium support.** Only Visium (FFPE + FF) and CosMx are wired in. Adding `spatialdata-io` parsers for those is a planned next step but would change the QC notebook materially.
- **Subtype-specific references.** Synovial sarcoma and angiosarcoma are well-covered by the bundled reference; rare subtypes (epithelioid sarcoma, alveolar soft-part sarcoma, etc.) inherit cell-state labels from the closest available subtype, which can over-smooth real heterogeneity.
- **Cross-sample batch correction.** Each sample is analyzed independently. A multi-sample Harmony or scVI integration is shown in the notebook but not part of the default Snakemake DAG.
- **Niche annotation.** Niche labels are unsupervised k-means on neighborhood composition. They are stable but not biologically named automatically.
- **Clinical metadata.** Survival, treatment history, and response annotations from the original studies are not harmonized into a common schema. If you need that, you'll want to write a small mapper per study.
- **CosMx cell-level analysis.** The CosMx code path treats cells the way Visium treats spots; we don't run cell-level segmentation refinement (`sopa`-based) by default.

If any of the above is a blocker for what you're trying to do, open an issue — most are scoped-out, not impossible.
