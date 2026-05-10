# Reference DBs for ARIBA

Each subdirectory here is a named DB built once with `ariba prepareref` and consumed by `ariba run`. The runner selects a DB via flag (e.g. `--kleb-virulence`).

## Layout

```
refs/
  <db_name>/
    inputs/           # vendored allele FASTAs (committed)
    metadata.tsv      # ARIBA prepareref metadata: gene, cluster, seq_type, var
    manifest.json     # source pkg + version, build date, file checksums
    prepareref_out/   # built artefact (gitignored)
```

## DBs

- `kleb_virulence/` — five Klebsiella virulence loci (ybt / clb / iuc / iro / rmp), vendored from `kleborate`. Built by `mag_rescue.pp.build_ariba_ref --kleb-virulence`.

## Build

ARIBA runs CD-HIT to cluster alleles as the first step of `prepareref`; we use its defaults.
