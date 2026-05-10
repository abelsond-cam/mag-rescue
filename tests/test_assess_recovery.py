"""Tests for tl/assess_recovery.

End-to-end CLI is exercised on HPC against real cohort outputs. Here we
cover the parsing + tally helpers with synthetic reports.
"""

from __future__ import annotations

from pathlib import Path

from mag_rescue.tl.assess_recovery import (
    DB_LOCI,
    DB_LOCI_NTOTAL,
    _free_text_labels,
    _per_sample_locus_genes,
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


def _row(free_text: str) -> str:
    """Build a synthetic ariba report row with the right column count + label."""
    cols = ["x"] * 31
    cols[-1] = free_text
    return "\t".join(cols) + "\n"


def _write_report(path: Path, free_texts: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(_ARIBA_HEADER)
        for ft in free_texts:
            fh.write(_row(ft))


def test_free_text_labels_dedupes(tmp_path: Path):
    p = tmp_path / "x.report.tsv"
    _write_report(p, ["ybt:ybtA", "ybt:ybtA", "iuc:iucA"])
    assert _free_text_labels(p) == {"ybt:ybtA", "iuc:iucA"}


def test_per_sample_locus_genes_classifies_by_regex(tmp_path: Path):
    rdir = tmp_path / "reports"
    _write_report(rdir / "S1.report.tsv", ["ybt:ybtA", "ybt:ybtE", "rmp:rmpA", "rmp:rmpA2"])
    _write_report(rdir / "S2.report.tsv", ["iro:iroB", "iro:iroC"])
    _write_report(rdir / "S3.report.tsv", [])  # no hits

    out = _per_sample_locus_genes(rdir, DB_LOCI["kleb_virulence"])
    assert set(out) == {"S1", "S2", "S3"}
    # S1: ybt has 2 genes, rmpADC has rmpA only, rmpa2 has rmpA2 only
    assert out["S1"]["ybt"] == {"ybt:ybtA", "ybt:ybtE"}
    assert out["S1"]["rmpADC"] == {"rmp:rmpA"}
    assert out["S1"]["rmpa2"] == {"rmp:rmpA2"}
    assert out["S1"]["iuc"] == set()
    # S2: only iro
    assert out["S2"]["iro"] == {"iro:iroB", "iro:iroC"}
    assert out["S2"]["ybt"] == set()
    # S3: nothing
    assert all(s == set() for s in out["S3"].values())


def test_assess_tally_with_no_bacotype(tmp_path: Path):
    run_dir = tmp_path / "CG39"
    rdir = run_dir / "reports"
    _write_report(rdir / "S1.report.tsv", ["ybt:ybtA"])
    _write_report(rdir / "S2.report.tsv", ["ybt:ybtA", "ybt:ybtE"])
    _write_report(rdir / "S3.report.tsv", [])

    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    assert n_total == 3
    by_name = {r.name: r for r in rows}
    # 2 of 3 samples have a ybt gene (S1, S2)
    assert by_name["ybt"].aribasr_n == 2
    assert by_name["ybt"].aribasr_pct == 2 / 3
    # neither has all 11
    assert by_name["ybt"].aribasr_full_n == 0
    # no bacotype data
    assert by_name["ybt"].bacotype_sr_pct is None


def test_render_markdown_runs_with_and_without_bacotype(tmp_path: Path):
    run_dir = tmp_path / "CG39"
    rdir = run_dir / "reports"
    _write_report(rdir / "S1.report.tsv", ["ybt:ybtA"])
    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    DB_LOCI_NTOTAL.clear()
    DB_LOCI_NTOTAL.update({l.name: l.n_total_genes for l in DB_LOCI["kleb_virulence"]})

    md = render_markdown(cohort="CG39", n_total=n_total, rows=rows)
    assert "CG39" in md
    assert "ARIBA SR" in md
    # No comparison columns when no bacotype data
    assert "Δ (ARIBA − Bacotype SR)" not in md

    # Now with bacotype: simulate a row by injecting via the tsv path
    bact = tmp_path / "complete_vs_sr_genomes" / "penetrance" / "CG39.csv"
    bact.parent.mkdir(parents=True)
    bact.write_text(
        "feature,p_val_corr,penetrance_ratio,complete_penetrance,sr_penetrance,locus_concordance,n_complete,n_sr,p_val\n"
        "Yersiniabactin_bsc,1.0,0.7,0.6,0.8,0.99,11,665,0.1\n"
    )
    n_total2, rows2 = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=bact)
    md2 = render_markdown(cohort="CG39", n_total=n_total2, rows=rows2)
    assert "Δ (ARIBA − Bacotype SR)" in md2


def test_render_tsv_columns(tmp_path: Path):
    run_dir = tmp_path / "X"
    rdir = run_dir / "reports"
    _write_report(rdir / "S1.report.tsv", ["iuc:iucA"])
    n_total, rows = assess(ariba_run_dir=run_dir, db="kleb_virulence", bacotype_penetrance_csv=None)
    DB_LOCI_NTOTAL.clear()
    DB_LOCI_NTOTAL.update({l.name: l.n_total_genes for l in DB_LOCI["kleb_virulence"]})
    tsv = render_tsv(cohort="X", n_total=n_total, rows=rows)
    header = tsv.splitlines()[0].split("\t")
    assert header[0] == "cohort"
    assert "ariba_sr_pct" in header
    assert "bacotype_sr_pct" in header
