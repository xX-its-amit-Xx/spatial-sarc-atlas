"""
Build a sarcoma-tuned scRNA reference by combining public sarcoma atlases.

Two ways to provide each reference (set in `config/references.yaml`):

  1. `accession: <GSE...>`
     The script downloads the GEO series' supplementary files and parses
     them into AnnData. It auto-detects two common scRNA formats:

       - 10x-style mtx triplet: `matrix.mtx(.gz)` + `barcodes.tsv(.gz)`
         + `(features|genes).tsv(.gz)`
       - Dense expression matrix: `*counts*.tsv(.gz)` / `*expression*.csv(.gz)`
         paired with a metadata file (`*metadata*` / `*annot*`)

     You tell the script which metadata column carries the cell-type label
     (`cell_type_column:`) and, when the series ships several metadata
     files, which one to pick (`metadata_pattern:` regex).

  2. `local_path: <path.h5ad>`
     Point at a pre-converted h5ad you built yourself. Bypasses the GEO
     downloader entirely. Useful for atlases hosted on Zenodo, CellxGene
     Census, or dbGaP where GEO isn't the canonical source.

Per-reference conversions are cached at
`data/references/cache/<id>.h5ad` so re-runs are fast.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

log = logging.getLogger("build-ref")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


# ---------------------------------------------------------------------------
# GEO download + extraction
# ---------------------------------------------------------------------------

def download_geo_supp(accession: str, out_dir: Path) -> Path:
    """Pull supplementary files for a GEO series and flatten/extract them."""
    try:
        import GEOparse
    except ImportError as exc:
        raise SystemExit("GEOparse is required to fetch GEO data") from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Fetching GEO series %s -> %s", accession, out_dir)
    gse = GEOparse.get_GEO(geo=accession, destdir=str(out_dir / "_meta"),
                           silent=True)
    gse.download_supplementary_files(directory=str(out_dir), download_sra=False)

    # GEOparse drops files under per-GSM subdirectories; flatten so the
    # reader sees a single directory.
    for nested in list(out_dir.glob("*/*")):
        if nested.is_file() and nested.parent != out_dir / "_meta":
            target = out_dir / nested.name
            if not target.exists():
                shutil.move(str(nested), str(target))

    # Extract any TAR archives in place (common for series-level bundles).
    for tar_path in list(out_dir.glob("*.tar")):
        log.info("Extracting %s", tar_path.name)
        with tarfile.open(tar_path) as tf:
            tf.extractall(out_dir, filter="data")

    return out_dir


# ---------------------------------------------------------------------------
# File-discovery helpers
# ---------------------------------------------------------------------------

def _find(files: list[Path], *patterns: str) -> Path | None:
    """Case-insensitive search; first pattern that hits, first file wins."""
    for pat in patterns:
        cre = re.compile(pat, re.IGNORECASE)
        for f in files:
            if cre.search(f.name):
                return f
    return None


def _compression(path: Path) -> str | None:
    return "gzip" if path.suffix == ".gz" else None


# ---------------------------------------------------------------------------
# scRNA matrix readers
# ---------------------------------------------------------------------------

def read_mtx_triplet(directory: Path):
    """Read a 10x-style mtx + barcodes + features triplet from a directory."""
    import scanpy as sc

    files = [f for f in directory.rglob("*") if f.is_file()]
    mtx = _find(files, r"matrix\.mtx(\.gz)?$")
    bc  = _find(files, r"barcodes\.tsv(\.gz)?$")
    ft  = _find(files, r"(features|genes)\.tsv(\.gz)?$")
    if not (mtx and bc and ft):
        return None

    log.info("Reading 10x-style triplet from %s", mtx.parent)
    adata = sc.read_mtx(mtx).T   # mtx is genes x cells; transpose to cells x genes
    adata.obs_names = pd.read_csv(bc, header=None, sep="\t",
                                  compression=_compression(bc))[0].astype(str).values
    feat = pd.read_csv(ft, header=None, sep="\t", compression=_compression(ft))
    gene_col = feat.columns[1] if feat.shape[1] > 1 else feat.columns[0]
    adata.var_names = feat[gene_col].astype(str).values
    adata.var_names_make_unique()
    return adata


def read_dense_matrix(path: Path):
    """Read a dense TSV/CSV expression matrix, autodetecting orientation."""
    import anndata as ad

    log.info("Reading dense matrix %s", path.name)
    sep = "\t" if ".tsv" in path.name.lower() else ","
    df = pd.read_csv(path, sep=sep, index_col=0, compression=_compression(path))

    # scRNA matrices have genes >> cells. If rows < cols, the file is genes×cells
    # and needs transposing to cells×genes.
    if df.shape[0] >= df.shape[1]:
        log.info("  shape %s — transposing to cells×genes", df.shape)
        df = df.T
    else:
        log.info("  shape %s — already cells×genes", df.shape)

    return ad.AnnData(
        X=df.values.astype(np.float32),
        obs=pd.DataFrame(index=df.index.astype(str)),
        var=pd.DataFrame(index=df.columns.astype(str)),
    )


def attach_metadata(adata, metadata_path: Path, cell_type_column: str):
    """Merge a cell-metadata table into adata.obs and set `cell_type`."""
    log.info("Attaching metadata %s", metadata_path.name)
    sep = "\t" if ".tsv" in metadata_path.name.lower() else ","
    meta = pd.read_csv(metadata_path, sep=sep, compression=_compression(metadata_path))

    # Heuristic: the index column is usually 'cell', 'Cell', 'barcode',
    # 'cell_id', or whatever the first column is named.
    candidates = ["cell", "Cell", "barcode", "Barcode", "cell_id", "CellID",
                  "NAME", meta.columns[0]]
    for idx_col in candidates:
        if idx_col in meta.columns:
            meta = meta.set_index(idx_col)
            break

    meta.index = meta.index.astype(str)
    shared = adata.obs_names.intersection(meta.index)
    if len(shared) == 0:
        raise SystemExit(
            f"No overlap between matrix cells and metadata index.\n"
            f"  matrix[:5] = {list(adata.obs_names[:5])}\n"
            f"  meta[:5]   = {list(meta.index[:5])}"
        )
    log.info("  matched %d/%d cells to metadata", len(shared), adata.n_obs)
    adata = adata[shared].copy()
    adata.obs = meta.loc[shared].copy()

    if cell_type_column not in adata.obs.columns:
        raise SystemExit(
            f"cell_type_column '{cell_type_column}' not found in metadata.\n"
            f"  Available columns: {list(adata.obs.columns)}"
        )
    adata.obs["cell_type"] = adata.obs[cell_type_column].astype(str)
    return adata


def read_geo_scrna(directory: Path, cell_type_column: str,
                   metadata_pattern: str | None = None,
                   counts_format: str = "auto"):
    """Auto-detect a scRNA dataset in `directory` and return an annotated AnnData."""
    files = [f for f in directory.rglob("*") if f.is_file()]

    # Pick the metadata file.
    if metadata_pattern:
        meta_path = _find(files, metadata_pattern)
    else:
        meta_path = _find(files, r"cell.?metadata", r"cell.?annot",
                          r"metadata", r"annotations", r"clinical")
    if not meta_path:
        raise SystemExit(
            f"No metadata file found in {directory}. "
            f"Set `metadata_pattern:` in references.yaml to specify one."
        )

    # Try the mtx triplet first unless caller forces dense.
    if counts_format in ("auto", "mtx10x"):
        adata = read_mtx_triplet(directory)
        if adata is not None:
            return attach_metadata(adata, meta_path, cell_type_column)
        if counts_format == "mtx10x":
            raise SystemExit(f"counts_format=mtx10x but no mtx triplet in {directory}")

    # Fall back to a dense matrix file. Exclude the metadata file from the search.
    candidates = [f for f in files if f != meta_path]
    count_path = _find(candidates, r"counts", r"expression", r"tpm",
                       r"\.tsv(\.gz)?$", r"\.csv(\.gz)?$")
    if count_path is None:
        raise SystemExit(f"No counts file found in {directory}.")

    adata = read_dense_matrix(count_path)
    return attach_metadata(adata, meta_path, cell_type_column)


# ---------------------------------------------------------------------------
# Label harmonization + per-type subsampling
# ---------------------------------------------------------------------------

def harmonize_labels(series: pd.Series, mapping: dict[str, list[str]]) -> pd.Series:
    reverse: dict[str, str] = {}
    for canonical, aliases in mapping.items():
        for alias in aliases:
            reverse[alias.lower()] = canonical
        reverse[canonical.lower()] = canonical
    return series.astype(str).str.lower().map(reverse).fillna("other")


def subsample_per_type(adata, label_col: str, n_per: int, seed: int):
    rng = np.random.default_rng(seed)
    keep: list = []
    for label, idx in adata.obs.groupby(label_col, observed=True).groups.items():
        idx = list(idx)
        if len(idx) <= n_per:
            keep.extend(idx)
        else:
            keep.extend(rng.choice(idx, size=n_per, replace=False))
    return adata[keep].copy()


# ---------------------------------------------------------------------------
# Per-reference orchestration
# ---------------------------------------------------------------------------

def load_one_reference(ref: dict, cache_dir: Path):
    """Return an AnnData for one reference, preferring local_path > cache > GEO."""
    import anndata as ad

    lp = ref.get("local_path")
    if lp and Path(lp).exists():
        log.info("Loading %s from local h5ad %s", ref["id"], lp)
        return ad.read_h5ad(lp)

    cached = cache_dir / f"{ref['id']}.h5ad"
    if cached.exists():
        log.info("Loading %s from cache %s", ref["id"], cached)
        return ad.read_h5ad(cached)

    accession = ref.get("accession")
    if not accession:
        raise SystemExit(
            f"Reference {ref['id']} has neither `accession` nor a usable "
            f"`local_path`. Add one or the other."
        )
    if "cell_type_column" not in ref:
        raise SystemExit(
            f"Reference {ref['id']} uses GEO download but is missing the "
            f"`cell_type_column:` field. Set it to the metadata column "
            f"that carries cell-type labels."
        )

    raw_dir = cache_dir / "raw" / ref["id"]
    download_geo_supp(accession, raw_dir)
    adata = read_geo_scrna(
        raw_dir,
        cell_type_column=ref["cell_type_column"],
        metadata_pattern=ref.get("metadata_pattern"),
        counts_format=ref.get("counts_format", "auto"),
    )
    cached.parent.mkdir(parents=True, exist_ok=True)
    adata.write(cached)
    log.info("Cached %s -> %s (%d cells × %d genes)",
             ref["id"], cached, adata.n_obs, adata.n_vars)
    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",         default="config/references.yaml", type=Path)
    parser.add_argument("--out-combined",   required=True, type=Path)
    parser.add_argument("--out-signatures", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path,
                        default=Path("data/references/cache"))
    args = parser.parse_args(argv)

    cfg = yaml.safe_load(open(args.config))
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    import anndata as ad
    import scanpy as sc

    refs = []
    for ref in cfg.get("references", []):
        try:
            a = load_one_reference(ref, args.cache_dir)
        except SystemExit as exc:
            log.warning("Skipping %s: %s", ref["id"], exc)
            continue
        except Exception as exc:
            log.warning("Skipping %s: unexpected %s: %s",
                        ref["id"], type(exc).__name__, exc)
            continue
        if "cell_type" not in a.obs:
            log.warning("Reference %s missing `cell_type` after load; skipping.",
                        ref["id"])
            continue
        a.obs["dataset"] = ref["id"]
        refs.append(a)

    if not refs:
        log.error("No usable references — see warnings above.")
        return 1

    combined = ad.concat(refs, join="inner", merge="same",
                         label="dataset", index_unique="-")
    combined.obs["cell_type_l1"] = harmonize_labels(
        combined.obs["cell_type"], cfg["label_harmonization"]
    ).astype("category")

    log.info("Combined reference: %d cells × %d genes, %d cell types",
             combined.n_obs, combined.n_vars,
             combined.obs["cell_type_l1"].nunique())

    args.out_combined.parent.mkdir(parents=True, exist_ok=True)
    combined.write(args.out_combined)

    sub_cfg = cfg.get("subsample", {})
    sub = subsample_per_type(combined, "cell_type_l1",
                             n_per=sub_cfg.get("cells_per_type", 2500),
                             seed=sub_cfg.get("random_state", 42))
    sc.pp.filter_genes(sub, min_cells=5)
    sub.write(args.out_signatures)
    log.info("Wrote signature input: %s", args.out_signatures)
    return 0


if __name__ == "__main__":
    sys.exit(main())
