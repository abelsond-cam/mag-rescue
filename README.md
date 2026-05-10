# mag-rescue

ARIBA-driven virulence/AMR profiling on Klebsiella short-read genomes.

The pipeline runs ARIBA per sample against a vendored reference DB, with the DB selected by flag (e.g. `--kleb-virulence`). Short reads are streamed in from ENA per sample, processed, and discarded — no bulk pre-staging.

The first DB is `kleb-virulence`: five Klebsiella virulence loci (ybt / clb / iuc / iro / rmp), seeded from [Kleborate](https://github.com/klebgenomics/Kleborate)'s reference alleles. ARIBA's first step (`ariba prepareref`) uses CD-HIT to cluster alleles — we use defaults.

## Quick start

```bash
pixi install              # solve env (osx-arm64 + linux-64 lock)
pixi run test             # smoke tests
pixi run lint             # ruff
```

For implementation guidance, see [CLAUDE.md](CLAUDE.md).

## Project layout

| Path | Purpose |
|------|---------|
| `src/mag_rescue/pp/` | Preprocessing: accession loading, ref DB build, fastq fetch |
| `src/mag_rescue/tl/` | Tools: parse ARIBA outputs, build summary tables |
| `src/mag_rescue/pl/` | Plotting: virulence cluster heatmaps |
| `refs/` | Vendored reference FASTAs and ARIBA metadata TSVs (one subdir per DB) |
| `slurm_scripts/` | Slurm array job wrappers for the HPC runner |
| `tests/` | pytest |

## Reference DBs

Each DB lives at `refs/<name>/` as a self-contained set of vendored inputs plus a built ARIBA artefact. Currently only `kleb_virulence` is built; future DBs (`amr`, `mlst`, …) plug in via the same flag pattern.

```
refs/kleb_virulence/
  inputs/<kleborate_module>/*.fasta   # vendored allele FASTAs (committed)
  inputs/<kleborate_module>/profiles.tsv  # MLST-style ST mappings (committed, future use)
  metadata.tsv                        # 6-col ARIBA prepareref input (committed)
  manifest.json                       # source pkg + version, build timestamp, sha256s (committed)
  prepareref_out/                     # ariba prepareref output, ~7.6 MB (gitignored — rebuildable)
```

`manifest.json` records the *installed* kleborate version (read via `importlib.metadata.version()`, so it tracks the conda package — not a stale `version.py` string), the build timestamp (UTC), and a SHA-256 for every vendored file. Diff two manifests to see what changed when refreshing from a new Kleborate release.

### ARIBA runs in an apptainer container

ARIBA itself is unmaintained and incompatible with modern pysam/samtools, so we run it inside the biocontainers/ariba:2.13.3 apptainer image (frozen pysam 0.15.0, samtools 1.9, etc.). One-off pull on HPC:

```bash
mkdir -p ~/rds/.../processed/mag_rescue/containers
cd ~/rds/.../processed/mag_rescue/containers
apptainer pull --name ariba_213.sif docker://quay.io/biocontainers/ariba:2.13.3--py36hfc679d8_0
```

See [CLAUDE.md](CLAUDE.md) "Env management" for the full rationale.

### Build / refresh the reference DB

Runs on HPC (Kleborate is linux-only):

```bash
ssh login.hpc.cam.ac.uk
cd ~/workspace/mag-rescue && git pull
pixi run -e refbuild python -m mag_rescue.pp.build_ariba_ref \
    --kleb-virulence \
    --ariba-sif ~/rds/.../processed/mag_rescue/containers/ariba_213.sif \
    --threads 4 [--force]
```

The script:
1. locates the pixi-installed `kleborate` package (in the `refbuild` env),
2. copies the 39 allele FASTAs (5 virulence loci × 6 modules — `rmp` is split across `klebsiella__rmst` and `klebsiella__rmpa2` upstream) plus their `profiles.tsv` files into `refs/kleb_virulence/inputs/`,
3. emits `metadata.tsv` (one row per allele, six tab-separated columns; ARIBA silently rejects rows with the wrong column count),
4. writes `manifest.json`,
5. runs `ariba prepareref` via apptainer with default CD-HIT settings (cdhit_min_id 0.9). Most of the 39 input genes collapse to one cluster, with a few outliers.

Then commit the refreshed `inputs/`, `metadata.tsv`, `manifest.json` from HPC.

### Refreshing when Kleborate updates

```bash
pixi update kleborate     # bumps the refbuild env's kleborate
git add pixi.lock && git commit -m "chore: bump kleborate" && git push

# on HPC, rebuild
git pull && pixi install -e refbuild
pixi run -e refbuild python -m mag_rescue.pp.build_ariba_ref \
    --kleb-virulence \
    --ariba-sif ~/rds/.../containers/ariba_213.sif \
    --force
git add refs/kleb_virulence/ && git commit -m "data: refresh kleb-virulence from kleborate vX.Y.Z"
git push
```

## License

GPL-3.0 — see [LICENSE](LICENSE).
