# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project purpose

mag-rescue runs [ARIBA](https://github.com/sanger-pathogens/ariba) on ~80,000 Klebsiella short-read genomes against curated allele dictionaries — recovering virulence (and later AMR/MLST) profiles directly from reads. The pipeline shape is parametrised by a `--<dbname>` flag; the first DB is `--kleb-virulence` (5 loci: ybt / clb / iuc / iro / rmp), seeded from [Kleborate](https://github.com/klebgenomics/Kleborate)'s reference alleles. ARIBA's `prepareref` step runs CD-HIT to cluster alleles — we use its defaults.

Sibling repo: [Bacotype](https://github.com/abelsond-cam/Bacotype) — pangenome / GPA / mobile-element analysis. mag-rescue consumes Bacotype's `metadata_final_curated_all_samples_and_columns.tsv` to source accessions: filter rows where `is_refseq == False` AND `kpsc_final_list == True`, then take the `run_accession` column. (`related_sr_accession` is reserved for an on-hold long-vs-short-read comparison and is **not** the column to use.)

## Commands

```bash
pixi install -e dev             # solve + install dev env (cross-platform lock)
pixi run -e dev test            # pytest
pixi run -e dev lint            # ruff check + format --check
pixi run -e dev fmt             # ruff format

# Build a reference DB (one-off; uses a separate `refbuild` env that has kleborate).
pixi run -e refbuild python -m mag_rescue.pp.build_ariba_ref \
    --kleb-virulence --ariba-sif <RDS>/.../containers/ariba_213.sif

# Extract a cohort accession list (HPC only — reads Bacotype's metadata).
pixi run -e dev python -m mag_rescue.pp.extract_accessions \
    --metadata <RDS>/.../final/metadata_final_curated_all_samples_and_columns.tsv \
    --outdir   <RDS>/.../processed/mag_rescue/kleb_virulence/<cohort>/accessions \
    --version v1 [--sublineage SL23 | --clonal-group CG39]

# Submit a Slurm array (uses parallel_ariba.py to compose sbatch).
pixi run -e dev python -m mag_rescue.pp.parallel_ariba submit \
    --db kleb_virulence --run-name <cohort> \
    --mag-rescue-root <RDS>/.../processed/mag_rescue \
    --repo-dir ~/workspace/mag-rescue \
    --ariba-sif <RDS>/.../containers/ariba_213.sif

# After completion: tally + compare to Bacotype's penetrance.
pixi run -e dev python -m mag_rescue.tl.assess_recovery \
    --ariba-run-dir <RDS>/.../mag_rescue/kleb_virulence/<cohort> \
    --cohort <cohort> \
    --bacotype-dir <RDS>/.../complete_vs_sr_genomes
```

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

Adding a new DB (`--amr`, `--mlst`, …): vendor FASTAs, write `metadata.tsv`, then register the DB name in two places:
- `pp/build_ariba_ref.py:DB_REGISTRY` (build-time module list)
- `tl/assess_recovery.py:DB_LOCI` (per-locus tally + Bacotype-feature mapping)

The runner (`pp/run_ariba.py`) is DB-agnostic — it just resolves `refs/<db>/prepareref_out/` from the `--db` flag.

## Pipeline shape

Four scripts wired by Slurm:

1. **`pp/extract_accessions.py`** (one-shot) — filter Bacotype's metadata to a 5-col TSV (`run_accession, r1_url, r2_url, r1_md5, r2_md5`). Optional `--sublineage` / `--clonal-group` for cohort runs.
2. **`pp/parallel_ariba.py submit`** — compose + sbatch the array against that list. Subcommands: `submit`, `status`, `retry` (re-submits failed indices via sacct).
3. **`slurm_scripts/ariba_array.sh`** — thin shim per array task; reads line N of the list, calls `pp/run_ariba.py`. icelake-himem partition, 4 cpus, 12 GB, 2 h.
4. **`pp/run_ariba.py`** — per-sample worker: `curl` R1+R2 from ENA → md5 verify → `apptainer exec ariba_213.sif ariba run` → atomic copy `report.tsv` (and optionally `assembled_seqs.fa.gz` + `assemblies.fa.gz` with `--detailed`) to RDS → cleanup scratch. Idempotent on `report.tsv` presence.
5. **`tl/assess_recovery.py`** — post-hoc: tally per-locus presence + gene completeness; if Bacotype's `complete_vs_sr_genomes/penetrance/<cohort>.csv` is around, emit comparison tables. Output: `<run-dir>/assessment/<cohort>_recovery.{md,tsv}`.

Outputs land at `<RDS>/processed/mag_rescue/<db>/<run-name>/{reports,sample_logs,assembled_seqs,assemblies,accessions,assessment}/`. Slurm stdout/stderr at `<RDS>/processed/mag_rescue/slurm_logs/`.

## Project state (2026-05-10)

- Phase 4 ref DB built: `refs/kleb_virulence/` (kleborate 3.2.4, 39 alleles → 55 CD-HIT clusters).
- CG39 cohort run complete: 612/665 samples (47 transient curl failures, 5 stale-md5 ENA mismatches). Initial findings: ARIBA rescues iro/rmpADC/rmpa2/clb (4–21× over Bacotype's existing SR detection); ybt and iuc match the existing rate. See `<RDS>/.../CG39/assessment/CG39_recovery.md`.
- SL23, CG307, CG340 arrays in flight (jobs 29161949 / 29161950 / 29161952).
- Hard rows (~2,000 mixed-platform / multi-run) deferred — see `<accessions>/kleb_short_reads_v1.skipped.tsv`.

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
