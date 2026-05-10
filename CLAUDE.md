# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project purpose

mag-rescue runs [ARIBA](https://github.com/sanger-pathogens/ariba) on ~80,000 Klebsiella short-read genomes against curated allele dictionaries — recovering virulence (and later AMR/MLST) profiles directly from reads. The pipeline shape is parametrised by a `--<dbname>` flag; the first DB is `--kleb-virulence` (5 loci: ybt / clb / iuc / iro / rmp), seeded from [Kleborate](https://github.com/klebgenomics/Kleborate)'s reference alleles. ARIBA's `prepareref` step runs CD-HIT to cluster alleles — we use its defaults.

Sibling repo: [Bacotype](https://github.com/abelsond-cam/Bacotype) — pangenome / GPA / mobile-element analysis. mag-rescue consumes Bacotype's `metadata_final_curated_all_samples_and_columns.tsv` to source accessions: filter rows where `is_refseq == False` AND `kpsc_final_list == True`, then take the `run_accession` column. (`related_sr_accession` is reserved for an on-hold long-vs-short-read comparison and is **not** the column to use.)

## Commands

```bash
pixi install                    # solve and install env (osx-arm64 + linux-64)
pixi run test                   # pytest
pixi run lint                   # ruff check + format --check
pixi run fmt                    # ruff format

pixi run python -m mag_rescue.pp.build_ariba_ref --kleb-virulence   # build ref DB (Phase 4, future)
```

Production scripts run on Slurm: edit knobs at the top of the relevant `slurm_scripts/*.sh`, then `sbatch`.

## HPC connection

- Host: `login.hpc.cam.ac.uk` (CSD3, user `dca36`). 8-hour SSH ControlMaster configured in `~/.ssh/config`.
- Code at `/home/dca36/workspace/mag-rescue` (siblings under `/home/dca36/workspace/`).
- Data under `/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david` — full storage map (four roots: `project_k`, `personal_rds`, `bacformer_rds`, `cold_storage`) in [`Bacotype/docs/data/hpc_storage_overview.md`](../Bacotype/docs/data/hpc_storage_overview.md). mag-rescue outputs land under `project_k/david/processed/mag_rescue/`.
- For code changes prefer `git commit` → `push` → `pull` on HPC over rsync.

## Package layout

scanpy-style modules:

| Module | Purpose |
|--------|---------|
| `src/mag_rescue/pp/` | Preprocessing — accession loading, reference DB build (`ariba prepareref`), fastq fetch (`prefetch`/`fasterq-dump`) |
| `src/mag_rescue/tl/` | Tools/analysis — parse ARIBA per-sample reports, collate into wide tables |
| `src/mag_rescue/pl/` | Plotting — virulence cluster heatmaps |

## Reference DBs

Each DB is a self-contained subdir under `refs/`:

```
refs/<db_name>/
  inputs/           # vendored allele FASTAs (committed)
  metadata.tsv      # ARIBA prepareref metadata: gene, cluster, seq_type, var
  manifest.json     # source pkg + version, build date, file checksums
  prepareref_out/   # built artefact (gitignored — rebuilt by `ariba prepareref`)
```

Adding a new DB (`--amr`, `--mlst`, …) means: vendor FASTAs, write `metadata.tsv`, register the flag in `pp/build_ariba_ref.py` and `pp/run_ariba.py`. The runner stays one shape.

## Pipeline shape

Per-sample worker (Phase 5):

1. Take an accession (SRR/ERR).
2. `prefetch` → `fasterq-dump --threads $SLURM_CPUS_PER_TASK` into `$SLURM_TMPDIR`.
3. `ariba run refs/<db>/prepareref_out/ <r1> <r2> $SLURM_TMPDIR/out`.
4. Copy report TSV to RDS.
5. Delete scratch fastqs.

Idempotent: skip if output TSV exists.

## Code style

- Line length: 120; numpy docstrings (enforced by `ruff pydocstyle`).
- Ruff: B, BLE, C4, D, E, F, I, RUF100, TID, UP, W (see `pyproject.toml` for ignores).
- Python 3.10–3.12 supported.

## Env management

Pixi is the single dependency manager. ARIBA's binary deps (bowtie2, samtools, mummer, cd-hit, spades) come from bioconda; `ariba` itself comes from PyPI. The lock file is multi-platform (osx-arm64 + linux-64) — same lock works locally and on HPC.

### Platform split — local mac is for editing, HPC runs the pipeline

`pixi.toml` declares two platforms: `osx-arm64` (local Apple Silicon dev) and `linux-64` (HPC + CI). The Python deps in the default `[dependencies]` table are cross-platform; the bioinformatics binaries and ARIBA itself are scoped to `[target.linux-64...]` only. Three reasons it has to be this way:

1. **ARIBA does not build on Apple Silicon.** ARIBA bundles C code that uses x86 SSE intrinsics (`__m64`, `mmintrin.h`). `pip install ariba` on osx-arm64 fails with clang errors like "invalid conversion between vector type `__m64` and integer type". So `ariba` lives under `[target.linux-64.pypi-dependencies]`.
2. **`kleborate` is linux-only on bioconda.** It depends on `stxtyper`, which currently has no osx-arm64 build. Since kleborate is only consumed at ref-build time (vendoring FASTAs into `refs/kleb_virulence/inputs/`), and that runs on HPC anyway, this is fine.
3. **`pymummer` (an ARIBA dep) needs `nucmer` on PATH at build time, not just install time.** Its `setup.py` probes for the MUMmer binaries before the wheel is built and aborts if they aren't found. With pip's default *build isolation*, the build venv is fresh and doesn't see the conda-installed `mummer` binaries — the build fails. We work around this with:

   ```toml
   [pypi-options]
   no-build-isolation = ["pymummer", "ariba"]
   ```

   Build isolation is then disabled for those two packages; pip uses the active pixi env (which has `mummer`) and the build sees `nucmer` on PATH. Without this entry, `pixi install -e dev` will fail on linux-64 with `Cannot install because some programs from the MUMer package not found`.
4. **ARIBA needs `setuptools` at runtime.** Its `__init__.py` does `from pkg_resources import get_distribution`. `pkg_resources` ships in `setuptools`, which modern Python conda envs don't bundle by default. So `setuptools` is an explicit `[target.linux-64.dependencies]` entry; without it `ariba version` fails with `ModuleNotFoundError: No module named 'pkg_resources'`.

Net effect: on macOS you can edit, lint, run unit tests, and import the package, but `pixi run` of anything that touches ARIBA only works on HPC. CI (linux-64) exercises the full pipeline.

### Bringing up the env on a fresh machine

```bash
# one-time pixi install (skip if already on PATH)
curl -fsSL https://pixi.sh/install.sh | bash

cd <repo>
pixi install -e dev   # solves and installs the dev env from pixi.lock
pixi run -e dev test
pixi run -e dev lint
```

On HPC, the same lockfile is used — pixi resolves the linux-64 entries from it. Pixi version on HPC must be ≥ the version that wrote the lockfile (lockfile schema v6 needs pixi 0.68+); update with `curl -fsSL https://pixi.sh/install.sh | bash` if older.
