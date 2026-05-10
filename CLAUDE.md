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

Pixi manages our **own** code's runtime (Python + pandas + biopython + tqdm + requests). **ARIBA itself runs in an apptainer (Singularity) container**, not in pixi. Likewise samtools, bowtie2, spades, mummer, cd-hit, pysam — they all live inside the container.

### Why a container for ARIBA?

ARIBA is unmaintained since ~2019 and is incompatible with the modern bioinformatics stack:

1. **ARIBA does not build on Apple Silicon.** ARIBA bundles C with x86 SSE intrinsics (`__m64`, `mmintrin.h`). Local mac dev cannot run the pipeline regardless — only HPC.
2. **`pymummer` (an ARIBA dep) probes for `nucmer` at build time.** pip's default build isolation hides conda binaries from the build venv, so `pip install ariba` fails.
3. **ARIBA's `__init__.py` imports `pkg_resources`** — removed from `setuptools` 80+ (Aug 2025).
4. **ARIBA's `samtools_variants.py` calls `pysam.mpileup('-t', '-u', '-v', ...)`.** Those options were removed in samtools 1.10 / pysam 0.16 (Jan-Apr 2020) when BCF/VCF output moved to bcftools. We can't pin pysam <0.16 in pixi because that version only has Python 2 wheels.

(1) is unavoidable. (2), (3), (4) compound: ARIBA only works against a frozen 2019-era env that pixi can't reach from Python 3.12.

The fix: **biocontainer `quay.io/biocontainers/ariba:2.13.3--py36hfc679d8_0`**, which ships:
- ariba 2.13.3
- pysam 0.15.0
- samtools 1.9 (bundled in pysam)
- bowtie2, spades, mummer, cd-hit, fermi-lite — all working versions

Pull once on HPC:

```bash
apptainer pull --name ariba_213.sif docker://quay.io/biocontainers/ariba:2.13.3--py36hfc679d8_0
```

The container is committed nowhere — it's ~200 MB and lives on RDS at `<RDS>/processed/mag_rescue/containers/ariba_213.sif`.

### pixi.toml layout

Three environments, two on the cross-platform `[dependencies]` table:

| env | features | purpose |
|---|---|---|
| `default` | (none) | runtime (our Python wrappers — pandas, biopython, requests, etc.) |
| `dev` | `dev` | + pytest, ruff, pre-commit |
| `refbuild` | `refbuild` (linux-64 only) | + kleborate, for the one-off ref-DB vendoring |

The `runtime` solve-group shared by default+dev keeps the dev env identical to production for the wrappers. The `refbuild` env is isolated so kleborate's transitive deps don't constrain runtime.

### Bringing up on a fresh machine

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # skip if already on PATH
cd <repo>
pixi install -e dev
pixi run -e dev test
pixi run -e dev lint
```

On HPC, the same lockfile applies (it's multi-platform: osx-arm64 + linux-64). Pixi version must be ≥0.68 (lockfile schema v6).

Container pull is a separate one-off:

```bash
mkdir -p ~/rds/.../processed/mag_rescue/containers
cd ~/rds/.../processed/mag_rescue/containers
apptainer pull --name ariba_213.sif docker://quay.io/biocontainers/ariba:2.13.3--py36hfc679d8_0
```

Then pass `--ariba-sif <path>` to `pp/run_ariba.py`, `pp/build_ariba_ref.py`, and `pp/parallel_ariba.py submit`.
