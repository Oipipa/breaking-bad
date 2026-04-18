# Volcano Plots - Run run_20260418_145455

## What was done

Volcano plots were generated from the full limma result tables (`*_all.csv`) using:

- `R/volcano_plots.R`
- Input run folder: `artifacts/deg_runs/run_20260418_145455`

The plotting script classifies probes as significant using the same strict rule as `R/diffex_genes.R`:

- adjusted p-value cutoff: `adj.P.Val < 0.05`
- fold-change cutoff: `|logFC| >= log2(1.2)` (about `0.263`)

For each contrast, the script writes both:

- `*.png` (quick viewing)
- `*.pdf` (publication-quality vector)

Summary table: `volcano_summary.csv`.

## What the results show in this run

From `volcano_summary.csv`:

- `Placebo_week12_vs_baseline`: 2 significant probes (1 up, 1 down)
- All other contrasts: 0 significant probes under the strict filter

This means the strict DEG calling is very conservative for this dataset/contrast setup. The volcano figures are still useful for visualizing trend structure even when strict hits are few.

## How to interpret the volcano figures

- X-axis: `logFC`
- Y-axis: `-log10(adj.P.Val)`
- Dashed vertical lines: fold-change threshold (`+/- log2(1.2)`)
- Dashed horizontal line: adjusted p-value threshold (`0.05`)
- Red: upregulated significant probes
- Blue: downregulated significant probes
- Grey: not significant by strict criteria

## Recommended usage for the project goal

- Keep strict `_sig` DEG files as high-confidence calls.
- Use these volcano plots (from `_all`) for reporting and visualization.
- Use ranked pathway enrichment results in `../pathways/` to capture biologically meaningful shifts when strict DEG counts are low.
