"""Pure-Python tests for the metadata + manifest assembly logic in build_ariba_ref.

The full build (vendoring kleborate, running ariba prepareref) is exercised on
HPC. Here we only cover the platform-independent helpers.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.pp.build_ariba_ref import (
    _build_metadata,
    _gene_basename,
    _read_fasta_headers,
)


def test_gene_basename_strips_kleborate_allele_suffix():
    assert _gene_basename("iucA_1") == "iucA"
    assert _gene_basename("rmpA2_42") == "rmpA2"
    assert _gene_basename("ybt_no_underscore") == "ybt_no"  # last underscore is the split point
    assert _gene_basename("nounder") == "nounder"


def test_read_fasta_headers_parses_simple_fasta(tmp_path: Path):
    fa = tmp_path / "iucA.fasta"
    fa.write_text(">iucA_1 some description\nACGT\n>iucA_2\nTGCA\n>iucA_3\nAAAA\n")
    assert _read_fasta_headers(fa) == ["iucA_1", "iucA_2", "iucA_3"]


def test_build_metadata_emits_one_row_per_sequence_with_cluster_label(tmp_path: Path):
    inputs = tmp_path / "inputs"
    (inputs / "klebsiella__abst").mkdir(parents=True)
    (inputs / "klebsiella__abst" / "iucA.fasta").write_text(">iucA_1\nACGT\n>iucA_2\nTGCA\n")
    (inputs / "klebsiella__rmpa2").mkdir(parents=True)
    (inputs / "klebsiella__rmpa2" / "rmpA2.fasta").write_text(">rmpA2_1\nACGT\n")

    tsv = _build_metadata(
        inputs,
        [("klebsiella__abst", "iuc"), ("klebsiella__rmpa2", "rmp")],
    )
    rows = [r.split("\t") for r in tsv.strip().splitlines()]
    assert rows == [
        ["iucA_1", "1", "0", ".", "iuc:iucA"],
        ["iucA_2", "1", "0", ".", "iuc:iucA"],
        ["rmpA2_1", "1", "0", ".", "rmp:rmpA2"],
    ]
