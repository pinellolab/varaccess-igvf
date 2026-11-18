# varaccess-igvf

Builds the released varACCESS accessibility tables for the IGVF submission and scores them
against the Safe-site negative-control panel.

The screen is varACCESS caQTL V3 in HepG2: 7,499 oligos of 181 nt over four replicates.

## Files

`varACCESS_per_element_accessibility.csv`, one row for each of the 7,428 oligos that pass QC.

| column | |
|---|---|
| `name` | oligo name from the library sheet |
| `category`, `subcategory` | library design annotation |
| `scored_positions` | editable positions kept by the context mask |
| `mean_depth` | mean read depth over those positions |
| `n_reps` | replicates contributing |
| `score_rep1..4` | per-replicate score |
| `score` | mean of the replicate scores |
| `sd_replicate` | standard deviation across replicates |
| `l2fc`, `pvalue`, `qvalue` | score referenced to the Safe panel, and its test |

`varACCESS_variant_accessibility_effects.csv`, one row per REF/ALT pair over every panel that
has allele pairs, 3,367 rows: 1,040 caQTL, 1,043 non-caQTL, 1,027 lipid-GWAS and 257 LDLR
promoter.

| column | |
|---|---|
| `pair_id` | key within this file |
| `panel` | `caQTL`, `noncaQTL`, `18loci` or `LDLRpro`, the per-element table's `category` |
| `variant_id`, `variant_type` | chrom:pos:ref:alt, and `snv`, `insertion` or `deletion` |
| `peak_id` | the published liver caPeak, caQTL and non-caQTL rows only |
| `test_family` | the set `qvalue` is corrected within |
| `oligo_ref`, `oligo_alt` | the two oligo names, as in the per-element table |
| `n_reps` | replicates with both alleles scored |
| `score_ref`, `score_alt` | replicate-mean score of each allele |
| `delta` | ALT minus REF, the replicate-mean difference |
| `se`, `t`, `df`, `ci_lo`, `ci_hi` | moderated standard error, statistic, degrees of freedom and 95% interval |
| `pvalue`, `qvalue` | two-sided P, and Benjamini-Hochberg q within `test_family` |
| `log2fc`, `log2fc_se` | mean over replicates of log2(ALT / REF), and its standard error |

The families are `caqtl_panels`, the caQTL and non-caQTL panels corrected together, and
`ldlr_snv`, `gwas_snv` and `gwas_indel`.

## Score

Editable positions are the TC and GA dinucleotide contexts of the 181 nt insert; positions below
depth 20 are dropped.
A context mask drops positions whose 5-mer reaches the 9 nt barcode, phrase positions 10 and 11,
in every oligo.
In both alleles of a caQTL or non-caQTL pair it also drops positions within 2 nt of the variant,
so the two alleles of a pair are scored on the same positions; the 18loci A1/A2 oligos and the
controls take the barcode rule only.

Per oligo per replicate the score is the unweighted mean of the per-position edit fraction over
the kept positions.
Per oligo it is the mean of the replicate scores.
Nothing is weighted, trimmed, shrunk or background-subtracted.

## Allelic effects

Within each replicate, ALT score minus REF score; `delta` is the mean of those differences.
`log2fc` is the same contrast on a ratio scale; the test is on `delta`, so `pvalue` and `qvalue`
rank with `delta` rather than with `log2fc`.
It is tested by a moderated t (Smyth 2004) on the per-replicate differences.
The replicate is the unit; sites and reads are never treated as independent observations.
A pair needs two replicates to be tested; the rest keep their `delta` and stay out of the
multiple-testing families.

The prior variance is fitted on the caQTL and non-caQTL pairs as splines in log pair depth, log
pair-mean score and log editable-position count, then frozen and evaluated at each lipid-GWAS and
LDLR pair's own covariates.
Those panels carry insertions and deletions as well as SNVs, so each pair's kept positions are
fixed per pair rather than per oligo, which makes the two alleles share the same genomic bases in
the same 5-mer context; for that reason `score_ref` and `score_alt` on those rows are not the same
quantity as the per-element table's `score` for the same oligo.

The 260 LDLR variants sit on 185 wild-type oligos, so variants sharing a wild-type oligo are
positively dependent; `oligo_ref` identifies the groups and the q values do not correct for it.

## Safe-site referencing

Each replicate is referenced to the median score of the Safe elements in that same replicate,
which cancels replicate differences in depth and global editing efficiency.

```
l2fc = mean over replicates of log2( score / Safe median score for that replicate )
```

The null is the same statistic on the Safe sites: centre is its median, sigma is a Gaussian scale
matched to its 5th/95th percentiles.

```
pvalue = 2 * norm.sf( |l2fc - centre| / sigma )
```

`qvalue` is Benjamini-Hochberg over the test elements only; Safe sites calibrate the null rather
than being hypotheses, so their `qvalue` is blank.
An element needs two replicates to be scored.

## Build

```sh
pip install -r requirements.txt
python prepare_release.py \
  --scores oligo_scores.csv.gz \
  --calls pair_calls.csv.gz --tests pair_tests.csv.gz \
  --annotation oligo_annotation.csv.gz --application application_pair_tests.csv.gz \
  --outdir released_data
```

The inputs are built by the varACCESS manuscript repository: `scripts/build_oligo_scores.py`,
`scripts/build_pairs.py` and `scripts/call_pairs.py` for the panels, and
`scripts/call_application_pairs.py` for the lipid-GWAS and LDLR pairs.
`accessibility_stats.py` is the Safe-site scoring on its own and can be run against a per-element
table that already has `name` and `score_rep<N>` columns:

```sh
python accessibility_stats.py in.csv out.csv
```
