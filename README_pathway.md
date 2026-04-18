# Ranked Metabolic Pathway Enrichment - Run run_20260418_145455

## What was done

Metabolic pathway enrichment was run with a ranked-list strategy (GSEA style), not ORA on strict DEG lists, using:

- `R/pathway_gsea_ranked.R`
- Input run folder: `artifacts/deg_runs/run_20260418_145455`

Method details:

1. Start from each contrast `*_all.csv` table.
2. Build a ranked gene list using limma moderated t-statistic (`t`), after collapsing probes to one entry per gene symbol (highest `|t|` kept).
3. Use metabolic pathway sets from MSigDB C2 (`msigdbr`) filtered by pathway-name keywords such as `METABOL`, `AMINO`, `LIPID`, `OXIDATIVE`, `TCA`, etc.
4. Run enrichment with `fgseaMultilevel`.
5. Mark significant pathways with `padj < 0.05`.

Outputs:

- `*_metabolic_gsea_all.csv` : all tested metabolic pathways
- `*_metabolic_gsea_sig.csv` : significant metabolic pathways (`padj < 0.05`)
- `metabolic_gsea_summary.csv` : per-contrast counts

## What the results show in this run

From `metabolic_gsea_summary.csv`:

- `Top_vs_Placebo_baseline`: 28 significant pathways
- `Top_week8_vs_baseline`: 28 significant pathways
- `Top_week12_vs_baseline`: 18 significant pathways
- `Interaction_week12`: 15 significant pathways
- `Interaction_week8`: 8 significant pathways
- `Top_vs_Placebo_week12`: 3 significant pathways
- `Top_vs_Placebo_week8`: 9 significant pathways
- `Placebo_week8_vs_baseline`: 7 significant pathways
- `Placebo_week12_vs_baseline`: 7 significant pathways

This is exactly why ranked enrichment is useful here: strict probe-level DEG files were sparse, but pathway-level signals are still detectable and interpretable.

## Examples of strongest signals

### Interaction week12 (`Interaction_week12_metabolic_gsea_sig.csv`)

Top significant pathways include:

- `REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY` (very strong positive NES)
- `REACTOME_SELENOAMINO_ACID_METABOLISM`
- `REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES`

### Top vs Placebo week12 (`Top_vs_Placebo_week12_metabolic_gsea_sig.csv`)

Top significant pathways include:

- `REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY`
- `REACTOME_SELENOAMINO_ACID_METABOLISM`
- `REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES`

## How to interpret columns

- `NES` : normalized enrichment score (direction and effect size)
- `padj` : multiple-testing adjusted p-value
- `leadingEdge` : core genes driving pathway enrichment

Interpretation rule of thumb:

- Higher `|NES|` and lower `padj` indicate stronger evidence.
- Positive NES means the pathway tends to be enriched among up-ranked genes.
- Negative NES means enrichment among down-ranked genes.

## Notes for paper-reproduction workflow

- Keep strict DEGs for high-confidence gene-level reporting.
- Use these ranked pathway results for biological interpretation and pathway sections.
- Prioritize interaction contrasts for treatment-effect-over-time interpretation.
