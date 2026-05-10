# Accession lists

This directory holds **only** the small audit summaries for accession-list extractions. The actual list TSVs (which can be ~MBs) live on RDS at `<project_k>/david/processed/mag_rescue/<db>/<run>/accessions/` and are not committed.

To regenerate:

```bash
ssh login.hpc.cam.ac.uk
cd ~/workspace/mag-rescue && git pull
pixi run -e dev python -m mag_rescue.pp.extract_accessions \
    --metadata ~/rds/rds-floto-bacterial-4k08a2yyQLw/david/final/metadata_final_curated_all_samples_and_columns.tsv \
    --outdir ~/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mag_rescue/kleb_virulence/all/accessions \
    --version v1
# then scp <RDS>/.../<run>/accessions/kleb_short_reads_v1.summary.txt back to this dir
```

The summary file records: row counts, skip reasons, and the source metadata path. Diffing two summaries tells you what changed when re-extracting.
