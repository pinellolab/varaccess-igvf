"""Build the two released tables from the varACCESS analysis tables."""

import argparse
from pathlib import Path

import pandas as pd

from accessibility_stats import score_table

ELEMENT_COLUMNS = {
    "sequence_name": "name",
    "category": "category",
    "subcategory": "subcategory",
    "n_pos": "scored_positions",
    "mean_depth": "mean_depth",
    "n_reps": "n_reps",
    "score_rep1": "score_rep1",
    "score_rep2": "score_rep2",
    "score_rep3": "score_rep3",
    "score_rep4": "score_rep4",
    "score": "score",
    "sd_replicate": "sd_replicate",
}

#: Panel names shared with the per-element table's ``category`` column.
PANEL_NAME = {"ldlr": "LDLRpro", "gwas": "18loci"}

APPLICATION_COLUMNS = {
    "variant_id": "variant_id",
    "variant_type": "variant_type",
    "family": "test_family",
    "oligo_ref": "oligo_ref",
    "oligo_alt": "oligo_alt",
    "n_reps": "n_reps",
    "ref_mean": "score_ref",
    "alt_mean": "score_alt",
    "delta": "delta",
    "se": "se",
    "t": "t",
    "df": "df",
    "ci_lo": "ci_lo",
    "ci_hi": "ci_hi",
    "p": "pvalue",
    "q": "qvalue",
    "log2fc": "log2fc",
    "log2fc_se": "log2fc_se",
}

CALL_COLUMNS = {
    "pair_id": "pair_id",
    "variant_id": "variant_id",
    "panel": "panel",
    "peak_id": "peak_id",
    "ref": "score_ref",
    "alt": "score_alt",
    "n_reps": "n_reps",
    "log2fc": "log2fc",
    "log2fc_se": "log2fc_se",
}

TEST_COLUMNS = {
    "pair_id": "pair_id",
    "a_trend_est": "delta",
    "a_trend_se": "se",
    "a_trend_stat": "t",
    "a_trend_df": "df",
    "a_trend_ci_lo": "ci_lo",
    "a_trend_ci_hi": "ci_hi",
    "a_trend_p": "pvalue",
    "a_trend_q": "qvalue",
}

EFFECT_ORDER = [
    "pair_id", "panel", "variant_id", "variant_type", "peak_id", "test_family",
    "oligo_ref", "oligo_alt", "n_reps", "score_ref", "score_alt",
    "delta", "se", "t", "df", "ci_lo", "ci_hi", "pvalue", "qvalue",
    "log2fc", "log2fc_se",
]


def per_element(scores_path):
    df = pd.read_csv(scores_path)[list(ELEMENT_COLUMNS)].rename(columns=ELEMENT_COLUMNS)
    return score_table(df.sort_values("name", ignore_index=True))


def benchmark_pairs(calls_path, tests_path, annotation_path):
    """The caQTL and non-caQTL panels, tested together against one BH family."""
    calls = pd.read_csv(calls_path)[list(CALL_COLUMNS)].rename(columns=CALL_COLUMNS)
    tests = pd.read_csv(tests_path)[list(TEST_COLUMNS)].rename(columns=TEST_COLUMNS)
    df = calls.merge(tests, on="pair_id", validate="one_to_one")

    ann = pd.read_csv(annotation_path)[["pair_id", "allele", "sequence_name"]]
    names = ann.pivot(index="pair_id", columns="allele", values="sequence_name")
    df = df.join(names.rename(columns={"REF": "oligo_ref", "ALT": "oligo_alt"}), on="pair_id")

    df["variant_type"] = "snv"
    df["test_family"] = "caqtl_panels"
    return df


def application_pairs(tests_path):
    """The LDLR promoter and lipid-GWAS pairs, each tested within its own BH family."""
    df = pd.read_csv(tests_path)
    df = df[df.n_reps > 0].copy()
    panel = df.application.map(PANEL_NAME)
    assert panel.notna().all(), "unknown application"
    out = df[list(APPLICATION_COLUMNS)].rename(columns=APPLICATION_COLUMNS)
    out["panel"] = panel.values
    out["peak_id"] = pd.NA
    out["pair_id"] = [f"{a}_{i:04d}" for a, i in
                      zip(df.application, df.groupby("application").cumcount() + 1)]
    return out


def variant_effects(calls_path, tests_path, annotation_path, application_path):
    df = pd.concat(
        [benchmark_pairs(calls_path, tests_path, annotation_path),
         application_pairs(application_path)],
        ignore_index=True,
    )
    return df[EFFECT_ORDER].sort_values(["panel", "pair_id"], ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", required=True, help="oligo_scores.csv.gz")
    ap.add_argument("--calls", required=True, help="pair_calls.csv.gz")
    ap.add_argument("--tests", required=True, help="pair_tests.csv.gz")
    ap.add_argument("--annotation", required=True, help="oligo_annotation.csv.gz")
    ap.add_argument("--application", required=True, help="application_pair_tests.csv.gz")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", default="varACCESS")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    elements, report = per_element(args.scores)
    elements_path = outdir / f"{args.prefix}_per_element_accessibility.csv"
    elements.to_csv(elements_path, index=False, float_format="%.6g")
    print(elements_path)
    print("\n".join(report))

    pairs = variant_effects(args.calls, args.tests, args.annotation, args.application)
    pairs_path = outdir / f"{args.prefix}_variant_accessibility_effects.csv"
    pairs.to_csv(pairs_path, index=False, float_format="%.6g")
    print(pairs_path)
    print(f"{len(pairs)} pairs, " + ", ".join(
        f"{n} {p}" for p, n in pairs.panel.value_counts().items()))
    print("  at FDR 5% within each family: " + ", ".join(
        f"{int((g.qvalue < 0.05).sum())}/{len(g)} {f}"
        for f, g in pairs.groupby("test_family")))
    print(f"  {int(pairs.pvalue.isna().sum())} untested")


if __name__ == "__main__":
    main()
