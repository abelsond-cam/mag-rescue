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

## License

GPL-3.0 — see [LICENSE](LICENSE).
