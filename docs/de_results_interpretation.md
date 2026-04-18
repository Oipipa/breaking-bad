# Differential Expression

## What the Pipeline Is Doing

The current public-data workflow is:

1. Build baseline-normalized expression deltas per subject:

- week 8 delta = week8 - baseline
- week 12 delta = week12 - baseline

2. Run limma differential expression per week:

- contrast = Topiramate vs Placebo
- outputs: full DE tables, overlap table, and ranked pathway input files

## Key Output Files

Gene-level:

- `artifacts/limma/limma_comparison_summary.csv`
- `artifacts/limma/de_public_week8_topiramate_vs_placebo_limma.csv`
- `artifacts/limma/de_public_week12_topiramate_vs_placebo_limma.csv`

## What the Current Results Are Showing

From the current run:

- Week 8: 1817 genes at nominal p<0.05, 0 genes at FDR<0.15.
- Week 12: 1957 genes at nominal p<0.05, 0 genes at FDR<0.15.

## How to Interpret This Pattern

1. Gene-level vs pathway-level sensitivity:

- It is possible to have weak gene-level FDR while still seeing pathway-level enrichment.
- This usually indicates distributed, coordinated effects across many genes in a pathway.

## What Results Should Generally Show (Expected Pattern)

For this public no-label setup, expected outputs are:

- Non-trivial nominal DE counts in both weeks.
- Limited gene-level FDR hits due to multiple-testing burden and modest effect sizes.
- More interpretable signal at pathway level (especially immune/inflammatory or metabolism/stress pathways if present).
- A mix of stable and time-dependent pathways across week 8 and week 12.

## Recommended Reporting Language

Use wording like:

- "In public-data treatment-effect contrasts, pathway-level enrichment provided clearer signal than single-gene FDR thresholds."
- "Week 12 showed a stronger pathway-level signature than week 8."
- "Findings should be interpreted as treatment-associated pathway trends rather than responder-specific biomarkers."
