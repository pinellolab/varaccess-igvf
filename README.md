# varaccess-igvf

Scores varACCESS elements against the Safe-site negative control panel for the IGVF
submission. Adds `l2fc`, `pvalue` and `qvalue` to a per-element accessibility table.

Two datasets, scored separately: **2024** (6 replicates, `SAFE1_Rep1..SAFE3_Rep6`) and
**2026** (4 replicates, `rep1..rep4`).

## Usage

```sh
pip install -r requirements.txt
python accessibility_stats.py in.csv out.csv
```

In-place is fine (`in.csv in.csv`), and re-running on an output is idempotent.

## Input

CSV with one row per element and these columns:

| column | |
|---|---|
| `name` | Safe sites must match `Safe<n>` or `Safe<n>_A1`; exactly 300 are required |
| `<rep>_edit_rate` | `sum(edits) / sum(total)` over the element's editable positions |
| `<rep>_mean_edits` | `sum(edits) / editable_positions` |

Replicates are discovered from the `_edit_rate` suffix, so either dataset works without
configuration. `overall_*` is ignored. Everything else is passed through untouched.

## Method

Reads per position is recovered as `mean_edits / edit_rate`. A replicate is used if that is
at least 10 and its edit rate is non-zero; an element is scored if at least 2 replicates
are used.

Each replicate is referenced to the median edit rate of the 300 Safe sites in that same
replicate, which cancels replicate differences in depth and global editing efficiency:

```
l2fc = mean over used replicates of log2( rate / Safe median rate for that replicate )
```

The null is the same statistic on the Safe sites: centre is its median, sigma is a Gaussian
scale matched to its 5th/95th percentiles. Safe residuals are near-Gaussian and their tail is
slightly lighter than normal, so this is mildly conservative.

```
pvalue = 2 * norm.sf( |l2fc - centre| / sigma )
```

`qvalue` is Benjamini-Hochberg over the test elements only; Safe sites calibrate the null
rather than being hypotheses, so their `qvalue` is blank.

The Gaussian fit replaces a straight empirical tail because 300 Safe sites floor an empirical
p-value at 1/301. It is calibrated against the Safe panel itself: 4.3% (2024) and 4.7% (2026)
of Safe sites reach p<0.05 against a nominal 5%.

Run the two datasets separately. Their Safe baselines differ about 3.5x
(2024 0.036-0.047, 2026 0.012-0.014), so a shared null would be wrong.

## Caveats

p-values far into the tail are extrapolated from a null estimated on 300 points; ranking is
meaningful there, the absolute value is not.

Safe sites are GC-poorer than test elements (median 0.370 vs 0.49-0.50) and Safe `l2fc` rises
with GC, so part of the signal for GC-rich elements is base composition.

The 2024 assay sits close to saturation (Safe baseline 0.036-0.047), leaving little headroom;
its `l2fc` values are compressed and the p-values are correspondingly low-powered.
