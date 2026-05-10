#!/bin/bash
#SBATCH --job-name=mag-ariba
#SBATCH --partition=icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#
# ariba_array.sh — thin Slurm shim for the per-sample ARIBA worker.
#
# Designed to be sbatch'd by `pp/parallel_ariba.py submit`, but also runnable
# directly for manual single-sample tests.
#
# Required env (set via sbatch --export=ALL,VAR=...):
#   LIST_FILE   absolute path to kleb_short_reads_*.tsv (5-col TSV with header)
#   DB          DB name under refs/<DB>/prepareref_out/, e.g. "kleb_virulence"
#   RUN_DIR     absolute path to the run output dir on RDS, e.g.
#               <RDS>/processed/mag_rescue/kleb_virulence/all
#   REPO_DIR    absolute path to the cloned mag-rescue repo (where pixi.toml lives)
# Optional:
#   SUBSET_METADATA   absolute path to a one-column file of accessions; matching
#                     accessions get --detailed
#
# This script: parses line $SLURM_ARRAY_TASK_ID + 1 of $LIST_FILE (skipping
# the header), invokes the Python worker with $SLURM_TMPDIR as scratch, and
# exits with the worker's exit code.

set -euo pipefail

: "${LIST_FILE:?LIST_FILE not set}"
: "${DB:?DB not set}"
: "${RUN_DIR:?RUN_DIR not set}"
: "${REPO_DIR:?REPO_DIR not set}"
: "${SUBSET_METADATA:=}"
: "${SLURM_ARRAY_TASK_ID:?must run as a Slurm array job}"
: "${SLURM_CPUS_PER_TASK:=4}"

# CSD3 doesn't always set $SLURM_TMPDIR — fall back to $TMPDIR, then /tmp.
# This matches the Bacotype convention.
SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORKDIR="${SCRATCH_BASE}/mag-ariba_${SLURM_JOB_ID:-local}_task${SLURM_ARRAY_TASK_ID}"
trap 'rm -rf "${WORKDIR}"' EXIT

# Read line N+1 (skip header).
LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))
ROW=$(sed -n "${LINE_NO}p" "${LIST_FILE}")
if [[ -z "${ROW}" ]]; then
    echo "ERROR: empty row at line ${LINE_NO} of ${LIST_FILE}" >&2
    exit 1
fi

IFS=$'\t' read -r ACC R1_URL R2_URL R1_MD5 R2_MD5 <<<"${ROW}"
if [[ -z "${ACC}" || -z "${R1_URL}" || -z "${R2_URL}" || -z "${R1_MD5}" || -z "${R2_MD5}" ]]; then
    echo "ERROR: malformed row at line ${LINE_NO}: ${ROW}" >&2
    exit 1
fi

DETAILED_FLAG=""
if [[ -n "${SUBSET_METADATA}" && -f "${SUBSET_METADATA}" ]]; then
    if grep -Fxq "${ACC}" "${SUBSET_METADATA}"; then
        DETAILED_FLAG="--detailed"
    fi
fi

echo "==== task=${SLURM_ARRAY_TASK_ID} acc=${ACC} detailed=${DETAILED_FLAG:-no} ===="

cd "${REPO_DIR}"
exec pixi run -e dev python -m mag_rescue.pp.run_ariba \
    --db "${DB}" \
    --accession "${ACC}" \
    --r1-url "${R1_URL}" --r2-url "${R2_URL}" \
    --r1-md5 "${R1_MD5}" --r2-md5 "${R2_MD5}" \
    --workdir "${WORKDIR}" \
    --run-dir "${RUN_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    ${DETAILED_FLAG}
