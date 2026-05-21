"""
Run cell2location on one Visium sample, given a pre-built sarcoma scRNA
reference. Writes the deconvolved AnnData and a per-cell-type spatial figure.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

log = logging.getLogger("c2l")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def fit_regression(ref):
    from cell2location.models import RegressionModel

    RegressionModel.setup_anndata(ref, labels_key="cell_type_l1",
                                  batch_key="dataset" if "dataset" in ref.obs else None)
    mod = RegressionModel(ref)
    mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002)
    ref = mod.export_posterior(ref, sample_kwargs={"num_samples": 1000,
                                                    "batch_size": 2500})
    means = ref.varm["means_per_cluster_mu_fg"].copy()
    means.columns = [c.replace("means_per_cluster_mu_fg_", "")
                     for c in means.columns]
    return means


def fit_spatial(adata, signatures):
    from cell2location.models import Cell2location

    shared = signatures.index.intersection(adata.var_names)
    adata = adata[:, shared].copy()
    signatures = signatures.loc[shared]

    Cell2location.setup_anndata(adata, batch_key=None)
    mod = Cell2location(
        adata,
        cell_state_df=signatures,
        N_cells_per_location=10,
        detection_alpha=20,
    )
    mod.train(max_epochs=15000, batch_size=None, train_size=1)
    adata = mod.export_posterior(adata, sample_kwargs={"num_samples": 1000,
                                                       "batch_size": adata.n_obs})
    return adata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path,
                        help="scRNA reference with raw counts + cell_type_l1")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    args = parser.parse_args(argv)

    log.info("Loading spatial sample %s", args.input)
    adata = sc.read_h5ad(args.input)

    log.info("Loading reference signatures %s", args.signatures)
    ref = sc.read_h5ad(args.signatures)

    log.info("Fitting cell2location regression on reference ...")
    sig_df = fit_regression(ref)

    log.info("Running cell2location on the spatial sample ...")
    adata = fit_spatial(adata, sig_df)

    # Convenience columns: dominant cell type per spot and per-type fractions
    prop = adata.obsm["q05_cell_abundance_w_sf"].copy()
    prop.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in prop.columns]
    norm = prop.div(prop.sum(axis=1), axis=0).fillna(0.0)
    for col in norm.columns:
        adata.obs[f"fraction_{col}"] = norm[col].values
    adata.obs["dominant_celltype"] = norm.idxmax(axis=1).astype("category")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write(args.output)

    # Quick multi-panel figure
    types = list(norm.columns)[:9]
    n = len(types)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    for ax, t in zip(np.atleast_1d(axes).ravel(), types):
        sc.pl.spatial(adata, color=f"fraction_{t}", ax=ax, show=False,
                      title=t, cmap="viridis", size=1.3)
    fig.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, bbox_inches="tight")
    log.info("Done — wrote %s and %s", args.output, args.figure)
    return 0


if __name__ == "__main__":
    sys.exit(main())
