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
