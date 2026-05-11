"""Tests for tl/assess_recovery.

End-to-end CLI is exercised on HPC against real cohort outputs. Here we
cover the parsing + tally + render helpers with synthetic reports.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.tl.assess_recovery import (
    DB_LOCI,
    ReportRow,
    _find_penetrance_csv,
    _per_sample_locus_rows,
    _read_report,
    assess,
    render_markdown,
    render_tsv,
)

_ARIBA_HEADER = (
    "#ariba_ref_name\tref_name\tgene\tvar_only\tflag\treads\tcluster\tref_len\t"
    "ref_base_assembled\tpc_ident\tctg\tctg_len\tctg_cov\tknown_var\tvar_type\t"
    "var_seq_type\tknown_var_change\thas_known_var\tref_ctg_change\tref_ctg_effect\t"
    "ref_start\tref_end\tref_nt\tctg_start\tctg_end\tctg_nt\tsmtls_total_depth\t"
    "smtls_nts\tsmtls_nts_depth\tvar_description\tfree_text\n"
)


def _row(
    *,
    free_text: str,
    reads: int = 100,
    ref_len: int = 1000,
    ref_base_assembled: int | None = None,
    pc_ident: float = 100.0,
    ctg_len: int = 1200,
    ctg_cov: float = 50.0,
) -> str:
    """Build a synthetic ariba report row with sensible numeric columns."""
    if ref_base_assembled is None:
        ref_base_assembled = ref_len
    cols = (
        [
            free_text,  # ariba_ref_name placeholder
            free_text,  # ref_name placeholder
            "1",
            "0",
            "27",
            str(reads),
            "cluster_x",
            str(ref_len),
            str(ref_base_assembled),
            str(pc_ident),
            "cluster_x.ctg.1",
            str(ctg_len),
            str(ctg_cov),
        ]
        + ["."] * 17
        + [free_text]
    )
    return "\t".join(cols) + "\n"


def _write_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(_ARIBA_HEADER)
        for line in lines:
            fh.write(line)


def test_read_report_parses_numeric_columns(tmp_path: Path):
    p = tmp_path / "x.report.tsv"
    _write_report(
        p,
        [
            _row(free_text="ybt:fyuA", reads=2416, ref_len=2022, ref_base_assembled=2022, ctg_cov=131.4),
            _row(free_text="ybt:irp1", reads=11450, ref_len=9492, ref_base_assembled=9492, ctg_cov=169.6),
        ],
    )
    rows = _read_report(p)
    assert len(rows) == 2
    assert rows[0].free_text == "ybt:fyuA"
    assert rows[0].ref_len == 2022
    assert rows[0].ref_base_assembled == 2022
    assert rows[0].ctg_cov == 131.4
    assert rows[1].reads == 11450


def test_per_sample_locus_rows_classifies_by_regex(tmp_path: Path):
    rdir = tmp_path / "reports"
    _write_report(
        rdir / "S1.report.tsv",
        [
            _row(free_text="ybt:ybtA"),
            _row(free_text="ybt:ybtE"),
            _row(free_text="rmp:rmpA"),
            _row(free_text="rmp:rmpA2"),
        ],
    )
    _write_report(rdir / "S2.report.tsv", [])  # no hits

    out = _per_sample_locus_rows(rdir, DB_LOCI["kleb_virulence"])
    assert set(out) == {"S1", "S2"}
    s1_ybt = {r.free_text for r in out["S1"]["ybt"]}
    assert s1_ybt == {"ybt:ybtA", "ybt:ybtE"}
    s1_rmpADC = {r.free_text for r in out["S1"]["rmpADC"]}
    s1_rmpa2 = {r.free_text for r in out["S1"]["rmpa2"]}
    assert s1_rmpADC == {"rmp:rmpA"}
    assert s1_rmpa2 == {"rmp:rmpA2"}
    assert out["S2"]["ybt"] == []


def test_assess_concordance_full_and_half(tmp_path: Path):
    """Full concordance = all 11 ybt genes; ≥50% = ≥6 (ceil(11/2))."""
    run_dir = tmp_path / "X"
    rdir = run_dir / "reports"
    # S1: 11 ybt genes → full
    _write_report(rdir / "S1.report.tsv", [_row(free_text=f"ybt:g{i}") for i in range(11)])
    # S2: 6 ybt genes → ≥50% but not full
    _write_report(rdir / "S2.report.tsv", [_row(free_text=f"ybt:g{i}") for i in range(6)])
    # S3: 5 ybt genes → below 50%
    _write_report(rdir / "S3.report.tsv", [_row(free_text=f"ybt:g{i}") for i in range(5)])
    # S4: no hits
    _write_report(rdir / "S4.report.tsv", [])

    # Note: each unique gene-label is what counts. The locus regex matches
    # `ybt:g0`..`ybt:g10` since the spec only requires `^ybt:`.
    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    by_name = {r.name: r for r in rows}
    ybt = by_name["ybt"]
    assert n_total == 4
    assert ybt.ariba_n == 3  # S1, S2, S3 are positive
    assert ybt.ariba_full_n == 1  # only S1
    assert ybt.ariba_half_or_more_n == 2  # S1 and S2 (S3 has 5 < ceil(11/2)=6)


def test_assess_recovery_quality_metrics(tmp_path: Path):
    """Recovery aggregates: kb recovered, % ref recovered, mean ctg_cov."""
    run_dir = tmp_path / "X"
    rdir = run_dir / "reports"
    # Single sample with one iuc gene: ref_len=2000, ref_base_assembled=1000, ctg_cov=50
    _write_report(
        rdir / "S1.report.tsv",
        [_row(free_text="iuc:iucA", ref_len=2000, ref_base_assembled=1000, ctg_cov=50.0)],
    )
    _, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    iuc = next(r for r in rows if r.name == "iuc")
    # 1000 / 1000 = 1.0 kb
    assert iuc.mean_kb_recovered == 1.0
    # 1000 / 2000 = 0.5
    assert iuc.mean_pct_ref_recovered == 0.5
    assert iuc.mean_ctg_cov == 50.0


def test_find_penetrance_csv_direct(tmp_path: Path):
    bact = tmp_path / "complete_vs_sr_genomes"
    (bact / "penetrance").mkdir(parents=True)
    direct = bact / "penetrance" / "CG39.csv"
    direct.write_text("feature\n")
    assert _find_penetrance_csv(bact, "CG39") == direct


def test_find_penetrance_csv_sl_to_cg_fallback(tmp_path: Path):
    bact = tmp_path / "complete_vs_sr_genomes"
    (bact / "penetrance").mkdir(parents=True)
    cg = bact / "penetrance" / "CG23.csv"
    cg.write_text("feature\n")
    assert _find_penetrance_csv(bact, "SL23") == cg


def test_find_penetrance_csv_no_match_returns_none(tmp_path: Path):
    bact = tmp_path / "complete_vs_sr_genomes"
    (bact / "penetrance").mkdir(parents=True)
    assert _find_penetrance_csv(bact, "CG99999") is None
    assert _find_penetrance_csv(bact, "SL99999") is None


def test_render_markdown_with_bacotype_uses_new_headers(tmp_path: Path):
    run_dir = tmp_path / "X"
    rdir = run_dir / "reports"
    _write_report(rdir / "S1.report.tsv", [_row(free_text="ybt:ybtA")])

    bact = tmp_path / "complete_vs_sr_genomes" / "penetrance" / "CG39.csv"
    bact.parent.mkdir(parents=True)
    bact.write_text(
        "feature,p_val_corr,penetrance_ratio,complete_penetrance,sr_penetrance,locus_concordance,n_complete,n_sr,p_val\n"
        "Yersiniabactin_bsc,1.0,0.7,0.6,0.8,0.99,11,665,0.1\n"
    )

    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=bact)
    md = render_markdown(cohort="CG39", n_total=n_total, rows=rows, penetrance_source=bact.name)
    # New headers
    assert "| Locus | Complete | sr-MAG | ARIBA | Δ |" in md
    assert "| Locus | n positive | % full concordance | % ≥50% concordance | mean kb recovered | mean depth (×) |" in md
    # Reflects penetrance source
    assert "CG39.csv" in md


def test_render_markdown_without_bacotype_omits_comparison_columns(tmp_path: Path):
    run_dir = tmp_path / "X"
    (run_dir / "reports").mkdir(parents=True)
    _write_report(run_dir / "reports" / "S1.report.tsv", [_row(free_text="ybt:ybtA")])
    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    md = render_markdown(cohort="X", n_total=n_total, rows=rows)
    assert "| Locus | ARIBA | n |" in md
    assert "sr-MAG" not in md
    assert "Complete" not in md


def test_render_tsv_has_new_columns(tmp_path: Path):
    run_dir = tmp_path / "X"
    (run_dir / "reports").mkdir(parents=True)
    _write_report(run_dir / "reports" / "S1.report.tsv", [_row(free_text="iuc:iucA")])
    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    tsv = render_tsv(cohort="X", n_total=n_total, rows=rows)
    header = tsv.splitlines()[0].split("\t")
    for col in (
        "ariba_pct",
        "ariba_full_frac",
        "ariba_half_or_more_frac",
        "mean_kb_recovered",
        "mean_pct_ref_recovered",
        "mean_ctg_cov",
        "srmag_pct",
        "complete_pct",
    ):
        assert col in header, f"missing {col}"
    # Old names should be gone
    assert "bacotype_sr_pct" not in header
    assert "ariba_sr_pct" not in header


def test_read_report_skips_malformed_numeric_rows(tmp_path: Path):
    p = tmp_path / "x.report.tsv"
    # Manually craft one bad row + one good row
    bad = "\t".join(["x", "x", "1", "0", "27", "not_a_number"] + ["."] * 25) + "\n"
    good = _row(free_text="ybt:fyuA")
    with p.open("w") as fh:
        fh.write(_ARIBA_HEADER)
        fh.write(bad)
        fh.write(good)
    rows = _read_report(p)
    assert len(rows) == 1
    assert rows[0].free_text == "ybt:fyuA"


def test_report_row_dataclass_exists():
    """Smoke: ReportRow is exposed and has expected fields."""
    r = ReportRow(
        free_text="ybt:fyuA",
        ref_len=100,
        ref_base_assembled=80,
        pc_ident=99.5,
        ctg_len=120,
        ctg_cov=10.0,
        reads=50,
    )
    assert r.free_text == "ybt:fyuA"
