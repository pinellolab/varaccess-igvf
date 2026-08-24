import argparse

import numpy as np
import pandas as pd
from scipy import stats

MIN_READS_PER_POSITION = 10
MIN_REPS = 2
SAFE_PATTERN = r"Safe\d+(_A1)?$"


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    reps = [c[: -len("_edit_rate")] for c in df.columns
            if c.endswith("_edit_rate") and not c.startswith("overall_")]

    rate = df[[f"{r}_edit_rate" for r in reps]].to_numpy(float)
    edits = df[[f"{r}_mean_edits" for r in reps]].to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = edits / rate
    rate = np.where((rate > 0) & (depth >= MIN_READS_PER_POSITION), rate, np.nan)

    safe = df["name"].str.match(SAFE_PATTERN).to_numpy()
    if safe.sum() != 300:
        raise SystemExit(f"expected 300 Safe elements, matched {safe.sum()}")

    ref = np.nanmedian(rate[safe], axis=0)
    per_rep = np.log2(rate / ref)
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

    df["l2fc"] = np.round(l2fc, 4)
    df["pvalue"] = pvalue
    df["qvalue"] = qvalue
    df.to_csv(args.output, index=False, float_format="%.6g")

    print(f"{len(df)} elements, {safe.sum()} Safe, {np.isfinite(l2fc).sum()} scored")
    print(f"reference edit rate per replicate: " + "  ".join(f"{r}={v:.5f}" for r, v in zip(reps, ref)))
    print(f"null centre={centre:+.4f} sigma={sigma:.4f}")
    print(f"Safe elements at p<0.05: {np.nanmean(pvalue[safe] < 0.05):.1%} (nominal 5%)")
    print(f"FDR 5%: {int(((qvalue < 0.05) & (l2fc > 0)).sum())} up, {int(((qvalue < 0.05) & (l2fc < 0)).sum())} down")


if __name__ == "__main__":
    main()
