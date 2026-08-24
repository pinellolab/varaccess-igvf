# varACCESS accessibility stats

Scores varACCESS elements against the Safe-site negative control panel.
Adds `l2fc`, `pvalue` and `qvalue` to a per-element accessibility table.

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

Replicates are discovered from the `_edit_rate` suffix, so any number of them works
(`rep1..rep4` for caQTL V3, `SAFE1_Rep1..SAFE3_Rep6` for HepG2). `overall_*` is ignored.
Everything else is passed through untouched.

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
p-value at 1/301. It is calibrated against the Safe panel itself: 4.7% (caQTL V3) and 4.3%
(HepG2) of Safe sites reach p<0.05 against a nominal 5%.

Datasets must be scored separately. Safe baselines differ about 3.5x between HepG2 and
caQTL V3, so a shared null would be wrong.
