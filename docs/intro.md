# spatial-sarc-atlas

A reproducible reanalysis pipeline for publicly available 10x Visium and CosMx spatial transcriptomics data on sarcoma samples.

This site renders the analysis notebooks with figures embedded so you can scroll through the results without running the pipeline yourself. For the full pipeline (download → QC → deconvolution → niche analysis), see the [GitHub repository](https://github.com/xX-its-amit-Xx/spatial-sarc-atlas).

## Navigating this book

- **Analysis walkthroughs**: notebooks 01–05, each self-contained and re-runnable.
- **Cookbook**: copy-paste recipes for common follow-up tasks once you have processed AnnData files.
- **Limitations**: honest list of what this pipeline does *not* do yet.

## Citation

If you use any code or figures from this work, please cite the original studies whose data is reanalyzed (listed in [config/samples.yaml](https://github.com/xX-its-amit-Xx/spatial-sarc-atlas/blob/main/config/samples.yaml)) and the underlying method packages: scanpy, squidpy, cell2location, Tangram, scvi-tools.
