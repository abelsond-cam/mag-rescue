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
    _retry_loop,
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


def test_array_spec_indices_collapses_consecutive_runs():
    """Slurm rejects too-long --array= strings — range-collapse to fit under 4KB."""
    # Consecutive runs collapsed to N-M
    assert _array_spec_indices([1, 2, 3, 4, 5], concurrency=10) == "1-5%10"
    # Mix of singletons, pairs, and longer runs
    assert _array_spec_indices([3, 5, 6, 8, 9, 10, 50], concurrency=5) == "3,5-6,8-10,50%5"
    # Two-element run still emits N-M (not N,M)
    assert _array_spec_indices([7, 8], concurrency=1) == "7-8%1"
    # De-dup + sort
    assert _array_spec_indices([5, 3, 5, 4], concurrency=2) == "3-5%2"


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


# ---------------------------------------------------------------------------
# Auto-retry loop tests (pure logic; subprocess interaction injected via mocks)
# ---------------------------------------------------------------------------


def _states_factory(scripts: list[list[tuple[int, str]]]):
    """Build a states_fn that returns successive scripted state lists per call."""
    calls = {"i": 0}

    def states_fn(_jobid: str) -> list[tuple[int, str]]:
        i = calls["i"]
        result = scripts[min(i, len(scripts) - 1)]
        calls["i"] += 1
        return result

    return states_fn


def test_retry_loop_clear_on_first_wave():
    """Initial wave has zero failures → status=clear, no retry submitted."""
    submitted: list[list[int]] = []

    def retry_submit(failed: list[int]) -> str:
        submitted.append(failed)
        return f"job-retry-{len(submitted)}"

    result = _retry_loop(
        initial_jobid="job-1",
        max_retries=3,
        wait_fn=lambda _: None,
        states_fn=_states_factory([[(1, "COMPLETED"), (2, "COMPLETED")]]),
        retry_submit_fn=retry_submit,
    )
    assert result["status"] == "clear"
    assert len(result["waves"]) == 1
    assert submitted == []


def test_retry_loop_clears_after_retries():
    """Wave 0: 3 failed. Wave 1: 1 failed. Wave 2: clear. Returns clear."""
    submitted: list[list[int]] = []

    def retry_submit(failed: list[int]) -> str:
        submitted.append(failed)
        return f"job-retry-{len(submitted)}"

    result = _retry_loop(
        initial_jobid="job-1",
        max_retries=3,
        wait_fn=lambda _: None,
        states_fn=_states_factory(
            [
                [(1, "COMPLETED"), (2, "FAILED"), (3, "FAILED"), (4, "FAILED"), (5, "COMPLETED")],
                [(2, "COMPLETED"), (3, "FAILED"), (4, "COMPLETED")],
                [(3, "COMPLETED")],
            ]
        ),
        retry_submit_fn=retry_submit,
    )
    assert result["status"] == "clear"
    assert len(result["waves"]) == 3
    assert submitted == [[2, 3, 4], [3]]  # two retry submissions, with shrinking failed sets


def test_retry_loop_stops_on_no_progress():
    """If failed count doesn't shrink between waves, declare no_progress."""
    submitted: list[list[int]] = []

    def retry_submit(failed: list[int]) -> str:
        submitted.append(failed)
        return f"job-retry-{len(submitted)}"

    result = _retry_loop(
        initial_jobid="job-1",
        max_retries=5,
        wait_fn=lambda _: None,
        states_fn=_states_factory(
            [
                [(1, "FAILED"), (2, "FAILED"), (3, "COMPLETED")],
                [(1, "FAILED"), (2, "FAILED"), (3, "COMPLETED")],  # same failures
            ]
        ),
        retry_submit_fn=retry_submit,
    )
    assert result["status"] == "no_progress"
    assert len(result["waves"]) == 2
    assert submitted == [[1, 2]]  # one retry attempt before declaring no progress


def test_retry_loop_respects_max_retries():
    """If failures keep shrinking but never reach zero, stop at max_retries."""
    submitted: list[list[int]] = []

    def retry_submit(failed: list[int]) -> str:
        submitted.append(failed)
        return f"job-retry-{len(submitted)}"

    # Wave 0: 5 failed, wave 1: 4, wave 2: 3, wave 3: 2 (max=2 allows up to 3 waves: 0,1,2 cycle)
    result = _retry_loop(
        initial_jobid="job-1",
        max_retries=2,
        wait_fn=lambda _: None,
        states_fn=_states_factory(
            [
                [(i, "FAILED") for i in range(5)],
                [(i, "FAILED") for i in range(4)],
                [(i, "FAILED") for i in range(3)],
            ]
        ),
        retry_submit_fn=retry_submit,
    )
    assert result["status"] == "max_retries"
    assert len(result["waves"]) == 3  # initial + 2 retries
