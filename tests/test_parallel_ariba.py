"""Tests for the orchestrator's pure-function helpers.

The actual sbatch / sacct interactions are exercised on HPC.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.pp.parallel_ariba import (
    _array_spec_indices,
    _array_spec_range,
    _build_sbatch_cmd,
    _count_data_rows,
    _list_path,
    _parse_jobid,
    _run_dir,
    _slurm_logs_dir,
    _write_list_for_test,
)


def test_run_and_list_path_layout(tmp_path: Path):
    root = tmp_path / "mag_rescue"
    rd = _run_dir(root, "kleb_virulence", "all")
    assert rd == root / "kleb_virulence" / "all"
    lp = _list_path(root, "kleb_virulence", "all", "kleb_short_reads_v1.tsv")
    assert lp == rd / "accessions" / "kleb_short_reads_v1.tsv"
    assert _slurm_logs_dir(root) == root / "slurm_logs"


def test_count_data_rows_excludes_header(tmp_path: Path):
    p = _write_list_for_test(tmp_path / "list.tsv", n_data_rows=42)
    assert _count_data_rows(p) == 42


def test_array_spec_range_caps_to_n_samples():
    assert _array_spec_range(n_total=100, n_samples=-1, concurrency=10) == "1-100%10"
    assert _array_spec_range(n_total=100, n_samples=10, concurrency=10) == "1-10%10"
    # n_samples larger than total -> clamped
    assert _array_spec_range(n_total=5, n_samples=100, concurrency=10) == "1-5%10"


def test_array_spec_indices_handles_explicit_list():
    assert _array_spec_indices([3, 17, 42], concurrency=20) == "3,17,42%20"


def test_parse_jobid_finds_id():
    assert _parse_jobid("Submitted batch job 29154321\n") == "29154321"
    assert _parse_jobid("error: nothing here") is None


def test_build_sbatch_cmd_passes_all_env_vars(tmp_path: Path):
    cmd = _build_sbatch_cmd(
        repo_dir=Path("/repo"),
        list_path=tmp_path / "list.tsv",
        db="kleb_virulence",
        run_dir=tmp_path / "rd",
        slurm_logs_dir=tmp_path / "slurm_logs",
        array_spec="1-10%5",
        ariba_sif=Path("/sif/ariba.sif"),
        subset_metadata=None,
    )
    assert cmd[0] == "sbatch"
    assert any(a.startswith("--array=1-10%5") for a in cmd)
    assert any(a.startswith("--output=") and a.endswith("/mag-ariba_%A_%a.out") for a in cmd)
    export = next(a for a in cmd if a.startswith("--export="))
    # Required env vars all present in --export
    for var in (
        "ALL",
        "LIST_FILE=",
        "DB=kleb_virulence",
        "RUN_DIR=",
        "REPO_DIR=/repo",
        "ARIBA_SIF=/sif/ariba.sif",
        "SUBSET_METADATA=",
    ):
        assert var in export, f"missing {var} in {export}"
    assert cmd[-1].endswith("ariba_array.sh")


def test_build_sbatch_cmd_includes_subset_metadata_path(tmp_path: Path):
    vip = tmp_path / "vip.txt"
    cmd = _build_sbatch_cmd(
        repo_dir=Path("/repo"),
        list_path=tmp_path / "list.tsv",
        db="kleb_virulence",
        run_dir=tmp_path / "rd",
        slurm_logs_dir=tmp_path / "slurm_logs",
        array_spec="1-1%1",
        ariba_sif=Path("/sif/ariba.sif"),
        subset_metadata=vip,
    )
    export = next(a for a in cmd if a.startswith("--export="))
    assert f"SUBSET_METADATA={vip}" in export
