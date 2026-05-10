"""Tests for run_ariba's pure-function helpers.

End-to-end (curl, ariba run) is exercised on HPC. Here we only cover
md5 verification and atomic copy.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.pp.run_ariba import _atomic_copy, _md5sum, _verify_md5


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
