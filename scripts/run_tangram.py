"""
Run Tangram on one Visium sample. Maps scRNA cells onto spatial spots and
produces per-spot cell-type abundance estimates.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

log = logging.getLogger("tangram")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--figure", required=True, type=Path)
    parser.add_argument("--n-markers", type=int, default=100,
                        help="Top marker genes per cell type to use as training genes")
    args = parser.parse_args(argv)

    import tangram as tg

    log.info("Loading inputs ...")
    spatial = sc.read_h5ad(args.input)
    sc_ref  = sc.read_h5ad(args.reference)

    if "cell_type_l1" not in sc_ref.obs:
        raise SystemExit("Reference must contain `cell_type_l1` in .obs")

    log.info("Computing markers per cell type ...")
    sc.tl.rank_genes_groups(sc_ref, "cell_type_l1", method="wilcoxon",
                            n_genes=args.n_markers)
    markers = (
        sc.get.rank_genes_groups_df(sc_ref, group=None)
          .groupby("group")
          .head(args.n_markers)["names"]
          .unique()
          .tolist()
    )

    shared = [g for g in markers if g in spatial.var_names]
    log.info("Using %d marker genes (%d shared with spatial)", len(markers), len(shared))

    tg.pp_adatas(sc_ref, spatial, genes=shared)
    log.info("Mapping cells onto spots ...")
    ad_map = tg.map_cells_to_space(
        adata_sc=sc_ref, adata_sp=spatial,
        device="cpu", mode="cells", density_prior="rna_count_based",
        num_epochs=500,
    )

    tg.project_cell_annotations(ad_map, spatial, annotation="cell_type_l1")

    # Tangram stores annotation matrix in obsm["tangram_ct_pred"]
    pred = spatial.obsm["tangram_ct_pred"].copy()
    pred = pred.div(pred.sum(axis=1), axis=0).fillna(0.0)
    for col in pred.columns:
        spatial.obs[f"fraction_{col}"] = pred[col].values
    spatial.obs["dominant_celltype"] = pred.idxmax(axis=1).astype("category")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    spatial.write(args.output)

    types = list(pred.columns)[:9]
    n = len(types); cols = 3; rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    for ax, t in zip(np.atleast_1d(axes).ravel(), types):
        sc.pl.spatial(spatial, color=f"fraction_{t}", ax=ax, show=False,
                      title=t, cmap="viridis", size=1.3)
    fig.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, bbox_inches="tight")

    log.info("Done — wrote %s and %s", args.output, args.figure)
    return 0


if __name__ == "__main__":
    sys.exit(main())
