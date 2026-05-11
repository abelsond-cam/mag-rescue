"""Tests for run_ariba's pure-function helpers.

End-to-end (curl, ariba run) is exercised on HPC. Here we only cover
md5 verification, atomic copy, and the --detailed destination path mapping.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.pp.run_ariba import (
    ARIBA_DETAILED_FILES,
    _atomic_copy,
    _detailed_dest,
    _md5sum,
    _verify_md5,
)


def test_md5sum_matches_known_value(tmp_path: Path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello world\n")
    assert _md5sum(p) == "6f5902ac237024bdd0c176cb93063dc4"


def test_verify_md5_pass_and_fail(tmp_path: Path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"hello world\n")
    assert _verify_md5(p, "6f5902ac237024bdd0c176cb93063dc4")
    assert not _verify_md5(p, "0" * 32)


def test_atomic_copy_creates_parents_and_no_tmp_left(tmp_path: Path):
    src = tmp_path / "src.bin"
    src.write_bytes(b"abc")
    dest = tmp_path / "deep" / "nested" / "dest.bin"
    _atomic_copy(src, dest)
    assert dest.read_bytes() == b"abc"
    assert not (dest.parent / "dest.bin.tmp").exists()


def test_detailed_dest_maps_single_dot_filenames(tmp_path: Path):
    """assembled_seqs.fa.gz → <run>/assembled_seqs/<acc>.fa.gz"""
    out = _detailed_dest(tmp_path, "DRR1", "assembled_seqs.fa.gz", "assembled_seqs")
    assert out == tmp_path / "assembled_seqs" / "DRR1.fa.gz"


def test_detailed_dest_maps_multi_dot_filenames(tmp_path: Path):
    """debug.report.tsv → <run>/debug_reports/<acc>.report.tsv"""
    out = _detailed_dest(tmp_path, "DRR1", "debug.report.tsv", "debug_reports")
    assert out == tmp_path / "debug_reports" / "DRR1.report.tsv"


def test_detailed_dest_handles_log_clusters_gz(tmp_path: Path):
    """log.clusters.gz → <run>/cluster_logs/<acc>.clusters.gz"""
    out = _detailed_dest(tmp_path, "DRR1", "log.clusters.gz", "cluster_logs")
    assert out == tmp_path / "cluster_logs" / "DRR1.clusters.gz"


def test_ariba_detailed_files_includes_all_seven_minus_version():
    """Sanity check: every ARIBA per-sample output we said we'd keep is registered."""
    expected_filenames = {
        "assembled_seqs.fa.gz",
        "assembled_genes.fa.gz",
        "assemblies.fa.gz",
        "debug.report.tsv",
        "log.clusters.gz",
    }
    actual = {fname for fname, _ in ARIBA_DETAILED_FILES}
    assert actual == expected_filenames
