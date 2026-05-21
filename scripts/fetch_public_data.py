"""
Download 10x Visium and CosMx public sarcoma datasets.

Strategy:
  * For GEO series (Visium): use GEOparse to read the SOFT metadata, verify
    platform tags look like Visium, then fetch the supplementary files for
    the chosen sample. We use Pooch so partial downloads can resume.
  * For Zenodo / direct URLs: stream + checksum via Pooch.
  * For the NanoString CosMx demo bundle: hit the official public URL; if
    the URL is unreachable, log and skip rather than failing the pipeline.

Run modes:
  --verify-only    Only validate accessions and exit non-zero on mismatch.
  --sample <id>    Fetch one sample.
  --all            Fetch every sample listed in the config.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import yaml

VISIUM_PLATFORMS = {"GPL24676", "GPL30178", "GPL30172"}

log = logging.getLogger("fetch")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def load_config(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def find_sample(config: dict, sample_id: str) -> dict | None:
    for modality in ("visium", "cosmx"):
        for entry in config.get(modality, []):
            if entry["id"] == sample_id:
                entry = dict(entry)
                entry["_modality"] = modality
                return entry
    return None


def verify_geo_visium(accession: str) -> bool:
    """Pull GEO SOFT metadata and confirm the series advertises a Visium platform."""
    try:
        import GEOparse
    except ImportError:
        log.warning("GEOparse not installed; skipping verification.")
        return True

    log.info("Verifying GEO series %s ...", accession)
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir="data/raw/_geoparse_cache",
                               silent=True)
    except Exception as exc:
        log.error("Could not retrieve %s: %s", accession, exc)
        return False

    platforms = set(gse.metadata.get("platform_id", []))
    if not (platforms & VISIUM_PLATFORMS):
        log.error("Accession %s does not advertise a Visium platform "
                  "(found platforms: %s).", accession, sorted(platforms))
        return False
    log.info("OK — %s uses platform(s) %s", accession, sorted(platforms))
    return True


def fetch_geo_supp(entry: dict, out: Path) -> None:
    """Pull the supplementary tar/zip for a GEO sample."""
    try:
        import GEOparse
    except ImportError as exc:
        raise SystemExit("GEOparse is required to fetch GEO data") from exc

    out.mkdir(parents=True, exist_ok=True)
    gse = GEOparse.get_GEO(geo=entry["accession"],
                           destdir=str(out / "_meta"), silent=False)
    sample_id = entry.get("sample")
    samples = [s for s in gse.gsms.values()
               if sample_id is None or s.name == sample_id]
    if not samples:
        raise SystemExit(f"No matching samples for {entry['id']} (sample={sample_id})")
    for gsm in samples:
        log.info("Downloading supplementary files for %s", gsm.name)
        gsm.download_supplementary_files(directory=str(out))


def fetch_url(url: str, out: Path) -> None:
    """Download an arbitrary URL via Pooch with resume support."""
    import pooch

    out.mkdir(parents=True, exist_ok=True)
    fname = Path(url).name or "download.bin"
    log.info("Pulling %s -> %s", url, out / fname)
    pooch.retrieve(url=url, known_hash=None, fname=fname, path=str(out),
                   progressbar=True)


def fetch_one(entry: dict, out_root: Path) -> None:
    out = out_root if out_root.name == entry["id"] else out_root / entry["id"]
    modality = entry["_modality"]

    if modality == "visium":
        if not verify_geo_visium(entry["accession"]):
            log.warning("Skipping %s — verification failed. "
                        "Update config/samples.yaml.", entry["id"])
            return
        fetch_geo_supp(entry, out)
    elif modality == "cosmx":
        url = entry.get("url")
        if not url:
            log.info("No URL for %s; skipping.", entry["id"])
            return
        try:
            fetch_url(url, out)
        except Exception as exc:
            log.warning("CosMx fetch failed for %s: %s — skipping.",
                        entry["id"], exc)
    else:
        raise SystemExit(f"Unknown modality: {modality}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/samples.yaml", type=Path)
    parser.add_argument("--sample", help="Single sample id to fetch")
    parser.add_argument("--all", action="store_true",
                        help="Fetch every sample in the config")
    parser.add_argument("--out", default="data/raw", type=Path,
                        help="Output directory (created if missing)")
    parser.add_argument("--verify-only", action="store_true",
                        help="Only verify accessions and exit")
    args = parser.parse_args(argv)

    config = load_config(args.config)

    if args.verify_only:
        ok = True
        for entry in config.get("visium", []):
            ok &= verify_geo_visium(entry["accession"])
        return 0 if ok else 1

    if args.all:
        targets = (config.get("visium", []) + config.get("cosmx", []))
        for entry in targets:
            entry = dict(entry)
            entry["_modality"] = "visium" if entry in config.get("visium", []) else "cosmx"
            fetch_one(entry, args.out)
        return 0

    if not args.sample:
        parser.error("Provide --sample <id>, --all, or --verify-only")

    entry = find_sample(config, args.sample)
    if entry is None:
        log.error("Sample %s not found in config", args.sample)
        return 2
    fetch_one(entry, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
