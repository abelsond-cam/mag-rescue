"""Tests for the URL/md5 parsing logic in extract_accessions.

The end-to-end filter (against a real Bacotype metadata.tsv) is exercised
on HPC. Here we only cover the platform-independent helpers.
"""

from __future__ import annotations

from mag_rescue.pp.extract_accessions import (
    _classify_row,
    _https_url,
    _pick_paired,
    _split_semicolon,
)


def test_https_url_adds_scheme():
    assert _https_url("ftp.sra.ebi.ac.uk/vol1/fastq/x") == "https://ftp.sra.ebi.ac.uk/vol1/fastq/x"
    assert _https_url("https://example.com/x") == "https://example.com/x"
    assert _https_url("ftp://example.com/x") == "ftp://example.com/x"


def test_split_semicolon_handles_empty_and_whitespace():
    assert _split_semicolon("") == []
    assert _split_semicolon("a;b;c") == ["a", "b", "c"]
    assert _split_semicolon(" a ; ; b ") == ["a", "b"]


def test_pick_paired_two_urls():
    urls = ["ftp.x/SRR1_1.fastq.gz", "ftp.x/SRR1_2.fastq.gz"]
    md5s = ["aaa", "bbb"]
    assert _pick_paired(urls, md5s) == (
        "https://ftp.x/SRR1_1.fastq.gz",
        "https://ftp.x/SRR1_2.fastq.gz",
        "aaa",
        "bbb",
    )


def test_pick_paired_three_urls_picks_pair_skips_orphan():
    urls = [
        "ftp.x/SRR1.fastq.gz",  # orphan
        "ftp.x/SRR1_1.fastq.gz",
        "ftp.x/SRR1_2.fastq.gz",
    ]
    md5s = ["orphan_md5", "r1_md5", "r2_md5"]
    assert _pick_paired(urls, md5s) == (
        "https://ftp.x/SRR1_1.fastq.gz",
        "https://ftp.x/SRR1_2.fastq.gz",
        "r1_md5",
        "r2_md5",
    )


def test_pick_paired_three_urls_without_pair_returns_none():
    # No _1/_2 endings — can't disambiguate.
    urls = ["ftp.x/a.fastq.gz", "ftp.x/b.fastq.gz", "ftp.x/c.fastq.gz"]
    md5s = ["a", "b", "c"]
    assert _pick_paired(urls, md5s) is None


def test_pick_paired_four_or_more_returns_none():
    urls = [f"ftp.x/SRR1_{i}.fastq.gz" for i in range(4)]
    assert _pick_paired(urls, ["m"] * 4) is None


def test_pick_paired_url_md5_length_mismatch_returns_none():
    assert _pick_paired(["a", "b"], ["m"]) is None


def _row(**overrides):
    """Build a metadata row with sensible Illumina-paired defaults."""
    base = {
        "run_accession": "SRR123",
        "metadata.runs.instrument.platform": "ILLUMINA",
        "fastq_ftp": "ftp.x/SRR123_1.fastq.gz;ftp.x/SRR123_2.fastq.gz",
        "fastq_md5": "aaa;bbb",
    }
    base.update(overrides)
    return base


def test_classify_row_clean_illumina_row_includes():
    out, reason = _classify_row(_row())
    assert reason is None
    assert out == ["SRR123", "https://ftp.x/SRR123_1.fastq.gz", "https://ftp.x/SRR123_2.fastq.gz", "aaa", "bbb"]


def test_classify_row_missing_accession():
    _, reason = _classify_row(_row(run_accession=""))
    assert reason == "no_run_accession"


def test_classify_row_missing_fastq_ftp():
    _, reason = _classify_row(_row(fastq_ftp=""))
    assert reason == "no_fastq_ftp"


def test_classify_row_pure_pacbio_skipped():
    _, reason = _classify_row(_row(**{"metadata.runs.instrument.platform": "PACBIO_SMRT"}))
    assert reason == "non_illumina"


def test_classify_row_multi_platform_skipped_even_with_illumina():
    _, reason = _classify_row(_row(**{"metadata.runs.instrument.platform": "ILLUMINA || OXFORD_NANOPORE"}))
    assert reason == "multi_platform_or_multi_run"


def test_classify_row_unsupported_url_count():
    _, reason = _classify_row(
        _row(
            fastq_ftp="ftp.x/a.fastq.gz",
            fastq_md5="aaa",
        )
    )
    assert reason == "unsupported_url_count_1"


def test_filter_with_sublineage_and_clonal_group(tmp_path):
    """End-to-end filter test against a tiny synthetic metadata TSV."""
    import csv

    from mag_rescue.pp.extract_accessions import _filter_to_kleb_short_reads

    rows = [
        # header
        [
            "kpsc_final_list",
            "is_refseq",
            "run_accession",
            "metadata.runs.instrument.platform",
            "fastq_ftp",
            "fastq_md5",
            "Sublineage",
            "Clonal group",
        ],
        # SL23 / CG23 — match
        ["True", "False", "SRR_A", "ILLUMINA", "ftp.x/A_1.fastq.gz;ftp.x/A_2.fastq.gz", "m1;m2", "SL23", "CG23"],
        # SL23 / CG39 — sublineage match but wrong CG
        ["True", "False", "SRR_B", "ILLUMINA", "ftp.x/B_1.fastq.gz;ftp.x/B_2.fastq.gz", "m1;m2", "SL23", "CG39"],
        # SL15 / CG39 — sublineage doesn't match
        ["True", "False", "SRR_C", "ILLUMINA", "ftp.x/C_1.fastq.gz;ftp.x/C_2.fastq.gz", "m1;m2", "SL15", "CG39"],
        # SL23 / CG23 but is_refseq — excluded by the kpsc/refseq filter
        ["True", "True", "SRR_D", "ILLUMINA", "ftp.x/D_1.fastq.gz;ftp.x/D_2.fastq.gz", "m1;m2", "SL23", "CG23"],
    ]
    p = tmp_path / "meta.tsv"
    with p.open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)

    inc_all, _ = _filter_to_kleb_short_reads(p)
    assert {r[0] for r in inc_all} == {"SRR_A", "SRR_B", "SRR_C"}

    inc_sl, _ = _filter_to_kleb_short_reads(p, sublineage="SL23")
    assert {r[0] for r in inc_sl} == {"SRR_A", "SRR_B"}

    inc_cg, _ = _filter_to_kleb_short_reads(p, clonal_group="CG39")
    assert {r[0] for r in inc_cg} == {"SRR_B", "SRR_C"}

    inc_both, _ = _filter_to_kleb_short_reads(p, sublineage="SL23", clonal_group="CG23")
    assert {r[0] for r in inc_both} == {"SRR_A"}
