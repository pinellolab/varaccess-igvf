"""Reference per-element accessibility to the Safe-site panel and test it."""

import argparse
import re

import numpy as np
import pandas as pd
from scipy import stats

MIN_REPS = 2
MIN_SAFE = 100
SAFE_PATTERN = r"Safe\d+(_A1)?$"
REP_PATTERN = re.compile(r"^score_rep(\d+)$")


def bh(p):
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return q
    order = np.argsort(p[ok])
    ranked = np.minimum.accumulate((p[ok][order] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def score_table(df):
    """Add l2fc, pvalue and qvalue to a per-element table. Returns the table and a report."""
    reps = sorted((c for c in df.columns if REP_PATTERN.match(c)),
                  key=lambda c: int(REP_PATTERN.match(c).group(1)))
    if not reps:
        raise SystemExit("no score_rep<N> columns found")

    score = df[reps].to_numpy(float)
    score = np.where(score > 0, score, np.nan)

    safe = df["name"].str.match(SAFE_PATTERN).to_numpy()
    if safe.sum() < MIN_SAFE:
        raise SystemExit(f"{safe.sum()} Safe elements, need at least {MIN_SAFE}")

    ref = np.nanmedian(score[safe], axis=0)
    per_rep = np.log2(score / ref)
    n_used = np.isfinite(per_rep).sum(axis=1)
    with np.errstate(invalid="ignore"):
        l2fc = np.where(n_used >= MIN_REPS, np.nanmean(per_rep, axis=1), np.nan)

    null = l2fc[safe]
    null = null[np.isfinite(null)]
    centre = np.median(null)
    spread = null - centre
    sigma = (np.quantile(spread, 0.95) - np.quantile(spread, 0.05)) / (2 * stats.norm.ppf(0.95))

    pvalue = 2 * stats.norm.sf(np.abs(l2fc - centre) / sigma)
    qvalue = np.full(len(df), np.nan)
    test = ~safe & np.isfinite(pvalue)
    qvalue[test] = bh(pvalue[test])

    df = df.copy()
    df["l2fc"] = np.round(l2fc, 4)
    df["pvalue"] = pvalue
    df["qvalue"] = qvalue

    report = [
        f"{len(df)} elements, {safe.sum()} Safe, {np.isfinite(l2fc).sum()} scored",
        "reference score per replicate: " + "  ".join(f"{r}={v:.5f}" for r, v in zip(reps, ref)),
        f"null centre={centre:+.4f} sigma={sigma:.4f}",
        f"Safe elements at p<0.05: {np.nanmean(pvalue[safe] < 0.05):.1%} (nominal 5%)",
        f"FDR 5%: {int(((qvalue < 0.05) & (l2fc > 0)).sum())} up, "
        f"{int(((qvalue < 0.05) & (l2fc < 0)).sum())} down",
    ]
    return df, report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    df, report = score_table(pd.read_csv(args.input))
    df.to_csv(args.output, index=False, float_format="%.6g")
    print("\n".join(report))


if __name__ == "__main__":
    main()
