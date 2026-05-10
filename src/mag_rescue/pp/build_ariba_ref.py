"""Build an ARIBA reference DB from vendored allele FASTAs.

Currently supports `--kleb-virulence` (Kleborate's 5+1 virulence cluster modules:
ybt, clb, iuc, iro, rmp, rmpa2 — cluster-labelled here as ybt/clb/iuc/iro/rmp,
with rmpA2 grouped under cluster ``rmp`` but kept separable in the metadata
description column as ``rmp:rmpA2``).

Future DBs (``--amr``, ``--mlst``) follow the same pattern: register them in
``DB_REGISTRY``, vendor source FASTAs, generate metadata.tsv, run prepareref.

Usage
-----
    pixi run -e dev python -m mag_rescue.pp.build_ariba_ref --kleb-virulence
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib
import importlib.metadata
import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger("build_ariba_ref")

# Each DB: list of (kleborate module, cluster label) pairs.
DB_REGISTRY: dict[str, dict] = {
    "kleb_virulence": {
        "source_pkg": "kleborate",
        "modules": [
            ("klebsiella__ybst", "ybt"),
            ("klebsiella__cbst", "clb"),
            ("klebsiella__abst", "iuc"),
            ("klebsiella__smst", "iro"),
            ("klebsiella__rmst", "rmp"),
            ("klebsiella__rmpa2", "rmp"),
        ],
    },
}

REPO_ROOT = Path(__file__).resolve().parents[3]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _source_pkg_version(pkg) -> str:
    """Return the installed package version.

    Prefers ``importlib.metadata.version()`` (reads the dist-info recorded by
    the package manager — conda/pip), then falls back to ``pkg.__version__``,
    then to a ``version.py`` ``__version__`` global. The metadata path matters
    because some packages (e.g. kleborate) ship a ``version.py`` whose string
    drifts from the actual installed version.
    """
    try:
        return importlib.metadata.version(pkg.__name__)
    except importlib.metadata.PackageNotFoundError:
        pass
    if hasattr(pkg, "__version__"):
        return str(pkg.__version__)
    version_py = Path(pkg.__file__).resolve().parent / "version.py"
    if version_py.exists():
        ns: dict = {}
        exec(version_py.read_text(), ns)
        return str(ns.get("__version__", "unknown"))
    return "unknown"


def _vendor_module(src_data: Path, dest_module: Path) -> dict[str, str]:
    """Copy ``*.fasta`` and ``profiles.tsv`` from src into dest. Return {filename: sha256}."""
    dest_module.mkdir(parents=True, exist_ok=True)
    sums: dict[str, str] = {}
    for src in sorted(src_data.iterdir()):
        if src.suffix == ".fasta" or src.name == "profiles.tsv":
            dest = dest_module / src.name
            shutil.copy2(src, dest)
            sums[src.name] = _sha256(dest)
    return sums


def _read_fasta_headers(fasta: Path) -> list[str]:
    """Return sequence names (no leading ``>``) in FASTA order."""
    names: list[str] = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            names.append(line[1:].split()[0])
    return names


def _gene_basename(seq_name: str) -> str:
    """Strip Kleborate's allele-number suffix: ``iucA_3`` → ``iucA``."""
    return seq_name.rsplit("_", 1)[0]


def _build_metadata(inputs_root: Path, modules: list[tuple[str, str]]) -> str:
    """Emit the ARIBA prepareref metadata TSV.

    Six tab-separated columns per row, as expected by
    ``ariba.sequence_metadata.SequenceMetadata``:

        name, gene_or_noncoding, variant_only, variant, variant_id, free_text

    Every row here: gene=1, variant_only=0, variant='.', variant_id='.',
    free_text='<cluster>:<gene>'. Five-column rows are silently rejected by
    ARIBA, leaving sequences orphaned in the FASTA.
    """
    rows: list[str] = []
    for mod, cluster in modules:
        for fasta in sorted((inputs_root / mod).glob("*.fasta")):
            for seq_name in _read_fasta_headers(fasta):
                rows.append("\t".join([seq_name, "1", "0", ".", ".", f"{cluster}:{_gene_basename(seq_name)}"]))
    return "\n".join(rows) + "\n"


def _check_apptainer_on_path() -> None:
    if shutil.which("apptainer") is None and shutil.which("singularity") is None:
        logger.error("Neither `apptainer` nor `singularity` on PATH. ariba runs in a container.")
        sys.exit(1)


def build(db_name: str, *, ariba_sif: Path, force: bool = False, threads: int = 1) -> None:
    """Vendor source FASTAs, write metadata + manifest, run ``ariba prepareref``."""
    spec = DB_REGISTRY[db_name]
    db_root = REPO_ROOT / "refs" / db_name
    inputs_root = db_root / "inputs"
    metadata_path = db_root / "metadata.tsv"
    manifest_path = db_root / "manifest.json"
    prepareref_out = db_root / "prepareref_out"

    _check_apptainer_on_path()
    if not ariba_sif.is_file():
        logger.error("apptainer SIF missing: %s", ariba_sif)
        sys.exit(2)

    if prepareref_out.exists() and not force:
        logger.error("%s already exists; pass --force to rebuild", prepareref_out)
        sys.exit(1)
    if force and prepareref_out.exists():
        shutil.rmtree(prepareref_out)

    src_pkg = importlib.import_module(spec["source_pkg"])
    src_pkg_dir = Path(src_pkg.__file__).resolve().parent
    src_pkg_version = _source_pkg_version(src_pkg)
    logger.info("Vendoring %s from %s (%s v%s)", db_name, src_pkg_dir, spec["source_pkg"], src_pkg_version)

    file_sums: dict[str, dict[str, str]] = {}
    for mod, _cluster in spec["modules"]:
        src_data = src_pkg_dir / "modules" / mod / "data"
        if not src_data.exists():
            logger.error("source module data missing: %s", src_data)
            sys.exit(1)
        dest = inputs_root / mod
        if dest.exists():
            shutil.rmtree(dest)
        sums = _vendor_module(src_data, dest)
        file_sums[mod] = sums
        logger.info("  %-30s  %d files", mod, len(sums))

    metadata_path.write_text(_build_metadata(inputs_root, spec["modules"]))
    n_seqs = sum(1 for _ in metadata_path.read_text().splitlines())
    logger.info("Wrote %s (%d sequences)", metadata_path.relative_to(REPO_ROOT), n_seqs)

    manifest = {
        "db_name": db_name,
        "built_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "source_pkg": spec["source_pkg"],
        "source_pkg_version": src_pkg_version,
        "source_pkg_path": str(src_pkg_dir),
        "modules": [{"name": m, "cluster": c} for m, c in spec["modules"]],
        "file_sha256": file_sums,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    logger.info("Wrote %s", manifest_path.relative_to(REPO_ROOT))

    fasta_args: list[str] = []
    for mod, _cluster in spec["modules"]:
        for fasta in sorted((inputs_root / mod).glob("*.fasta")):
            fasta_args.extend(["-f", str(fasta)])
    inner_cmd = [
        "ariba",
        "prepareref",
        *fasta_args,
        "-m",
        str(metadata_path),
        "--threads",
        str(threads),
        str(prepareref_out),
    ]
    cmd = ["apptainer", "exec", "-B", f"{db_root}:{db_root}", str(ariba_sif), *inner_cmd]
    logger.info("Running ariba prepareref (containerised) → %s", prepareref_out.relative_to(REPO_ROOT))
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        logger.error("ariba prepareref failed (returncode %d)", res.returncode)
        sys.exit(res.returncode)
    logger.info("Build complete: %s", db_root.relative_to(REPO_ROOT))


def main() -> None:
    """CLI entrypoint."""
    ap = argparse.ArgumentParser(description="Build an ARIBA reference DB from vendored allele FASTAs.")
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument(
        "--kleb-virulence",
        action="store_const",
        dest="db",
        const="kleb_virulence",
        help="Build the Kleborate-derived virulence DB (5 loci: ybt/clb/iuc/iro/rmp).",
    )
    ap.add_argument("--ariba-sif", type=Path, required=True, help="Path to the ariba apptainer container.")
    ap.add_argument("--force", action="store_true", help="Rebuild even if prepareref_out/ exists.")
    ap.add_argument("--threads", type=int, default=1, help="Threads for cd-hit inside prepareref.")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    build(args.db, ariba_sif=args.ariba_sif, force=args.force, threads=args.threads)


if __name__ == "__main__":
    main()
