# Pathway Analysis Results Interpretation

## Purpose of This Document

This document explains how to interpret the pathway-analysis outputs generated from the GSE107015 Topiramate-vs-placebo workflow.

The pathway-analysis step uses the limma outputs as input:

- preranked GSEA uses the full ranked gene lists: `artifacts/limma/*_ranked.rnk`
- over-representation analysis uses filtered gene lists: `artifacts/limma/*_pathway_input_p05.csv` and `artifacts/limma/*_pathway_input_fdr015.csv`

The main pathway outputs are stored in:

- `artifacts/pathway/gsea_all_results.csv`
- `artifacts/pathway/ora_all_results.csv`
- `artifacts/pathway/pathway_analysis_summary.csv`

## Analysis Question

The main biological question is:

> Which biological pathways show coordinated expression changes in Topiramate-treated participants compared with placebo after baseline correction?

The contrasts are:

- week 8: `Topiramate - Placebo` on `(week8 - baseline)`
- week 12: `Topiramate - Placebo` on `(week12 - baseline)`

Because the expression data come from whole blood, pathway-level findings should be interpreted as treatment-associated blood transcriptional signatures, not as direct measurements of pathway activity in brain or other tissues.

## Primary Evidence: GSEA Results

Use `gsea_all_results.csv` as the main interpretation file.

Important columns:

| Column | Meaning |
|---|---|
| `contrast` | Which comparison was tested, usually week 8 or week 12. |
| `collection` | MSigDB collection used, such as Hallmark, Reactome, or KEGG. |
| `pathway` | Name of the enriched gene set. |
| `pathway_description` | Short pathway description from MSigDB, when available. |
| `NES` | Normalized enrichment score. This gives direction and relative strength. |
| `padj` | FDR-adjusted pathway p-value. This is the main significance column. |
| `size` | Number of genes from the pathway included in the tested ranked list. |
| `leading_edge` | Genes driving the enrichment signal. |

### Direction of Enrichment

In this project:

- `NES > 0` means the pathway is enriched among genes with higher Topiramate-vs-placebo baseline-normalized expression change.
- `NES < 0` means the pathway is enriched among genes with lower Topiramate-vs-placebo baseline-normalized expression change.

Suggested wording:

- Positive NES: "The pathway was positively enriched in the Topiramate-vs-placebo contrast."
- Negative NES: "The pathway was negatively enriched in the Topiramate-vs-placebo contrast."

Avoid saying that a pathway is definitely "activated" or "inhibited" unless the pathway biology and leading-edge genes strongly support that interpretation.

## Secondary Evidence: ORA Results

Use `ora_all_results.csv` as supporting evidence.

Important columns:

| Column | Meaning |
|---|---|
| `contrast` | Week 8 or week 12 comparison. |
| `collection` | MSigDB collection used. |
| `threshold` | Which gene list was used, usually `p05` or `fdr015`. |
| `pathway` | Pathway/gene set name. |
| `overlap_size` | Number of significant genes overlapping the pathway. |
| `significant_size` | Number of genes in the tested foreground list. |
| `pathway_size` | Number of pathway genes in the test universe. |
| `gene_ratio` | Proportion of foreground genes overlapping the pathway. |
| `background_ratio` | Proportion of background genes in the pathway. |
| `padj` | FDR-adjusted ORA p-value. |
| `overlap_genes` | Genes responsible for the ORA overlap. |

ORA is more sensitive to the chosen gene cutoff than GSEA. Treat ORA as confirmatory or supportive rather than the main result.

## Recommended Significance Thresholds

Use this interpretation scheme:

| Result | Interpretation |
|---|---|
| `padj < 0.05` | Strong pathway-level evidence. |
| `padj < 0.25` | Exploratory pathway-level evidence, useful for hypothesis generation. |
| `padj >= 0.25` | Not enough evidence to report as an enriched pathway. |

Use the stricter `padj < 0.05` threshold for final claims whenever possible. If no pathways pass 0.05, report `padj < 0.25` findings clearly as exploratory.

## Step-by-Step Interpretation Workflow

### 1. Inspect the pathway summary

Open:

```text
artifacts/pathway/pathway_analysis_summary.csv
```

Check:

- whether any GSEA results have `n_fdr_lt_0_05 > 0`
- whether any GSEA results have `n_fdr_lt_0_25 > 0`
- whether results differ between week 8 and week 12
- whether ORA supports the same pathway categories

### 2. Identify top GSEA pathways

Open:

```text
artifacts/pathway/gsea_all_results.csv
```

Filter to:

```text
padj < 0.25
```

Then sort by:

```text
contrast, padj, abs(NES)
```

Create a short table with the top pathways:

| Contrast | Collection | Pathway | NES | padj | Direction | Main leading-edge genes |
|---|---|---:|---:|---:|---|---|
| week 8 |  |  |  |  |  |  |
| week 12 |  |  |  |  |  |  |

### 3. Group similar pathways into biological themes

Do not interpret every pathway as an independent discovery. Many gene sets overlap.

Group related pathways into themes such as:

- immune or inflammatory signalling
- cytokine or interferon response
- oxidative stress or reactive oxygen species biology
- mitochondrial or metabolic processes
- cell cycle or proliferation
- apoptosis or stress-response signalling
- signal transduction pathways

Example interpretation pattern:

```text
Several enriched pathways pointed to a shared immune/inflammatory theme rather than isolated independent findings. This was supported by enrichment of [PATHWAY 1], [PATHWAY 2], and [PATHWAY 3], with leading-edge genes including [GENE LIST].
```

### 4. Check whether week 8 and week 12 agree

Use the `contrast` column to compare timepoints.

Possible patterns:

| Pattern | Interpretation |
|---|---|
| Same pathway/theme in week 8 and week 12, same NES direction | More consistent treatment-associated signature. |
| Pathway significant only at week 12 | Possible delayed transcriptional response or stronger late effect. |
| Pathway significant only at week 8 | Possible early or transient response. |
| Opposite NES direction between weeks | Treat cautiously; may reflect time-dependent biology, noise, or unstable signal. |
| No pathways below FDR threshold | Limited evidence for robust pathway-level changes under this analysis design. |

### 5. Use ORA as supporting evidence

Open:

```text
artifacts/pathway/ora_all_results.csv
```

Filter to:

```text
padj < 0.05
```

Then ask:

- Do ORA pathways match the same biological themes as GSEA?
- Are the same genes appearing in `overlap_genes` and `leading_edge`?
- Are ORA findings only present for the nominal `p05` list, or also for the stricter `fdr015` list?

Stronger evidence:

```text
A pathway theme appears in both GSEA and ORA, with overlapping driver genes and the same timepoint pattern.
```

Weaker evidence:

```text
A pathway appears only in ORA from the nominal p<0.05 gene list and is not supported by GSEA.
```

## Results Summary Template

Fill this section after reviewing the CSV files.

### GSEA Summary

```text
At week 8, [NUMBER] pathways passed FDR < 0.05 and [NUMBER] pathways passed FDR < 0.25. The strongest positively enriched pathways were [PATHWAYS], suggesting higher Topiramate-vs-placebo baseline-normalized expression changes in genes related to [BIOLOGICAL THEME]. The strongest negatively enriched pathways were [PATHWAYS], suggesting lower Topiramate-vs-placebo expression changes in genes related to [BIOLOGICAL THEME].

At week 12, [NUMBER] pathways passed FDR < 0.05 and [NUMBER] pathways passed FDR < 0.25. The strongest positively enriched pathways were [PATHWAYS], while the strongest negatively enriched pathways were [PATHWAYS]. Compared with week 8, the week 12 results were [stronger/weaker/similar/different], indicating [INTERPRETATION].
```

### ORA Summary

```text
ORA provided [supporting / limited / no] evidence for the GSEA findings. The strongest ORA-supported themes were [BIOLOGICAL THEMES], based on enrichment of [PATHWAYS] among genes passing [p05/fdr015] thresholds. These results were driven by overlapping genes including [GENES].
```

### Integrated Interpretation

```text
Overall, the pathway results suggest that Topiramate treatment was associated with coordinated transcriptional changes in [BIOLOGICAL THEME 1], [BIOLOGICAL THEME 2], and [BIOLOGICAL THEME 3] in whole blood. The most consistent evidence came from [TIMEPOINT] and was supported by [GSEA only / GSEA and ORA]. Because these analyses are based on blood gene-expression changes and pathway-level enrichment, the findings should be interpreted as treatment-associated molecular signatures rather than direct proof that Topiramate activates or inhibits these pathways.
```

## If No Pathways Are Significant

Use this wording if no pathways pass the selected FDR threshold:

```text
No pathways passed the selected FDR threshold in the Topiramate-vs-placebo pathway analysis. This suggests that, under the current public-data analysis design, there is limited evidence for robust coordinated pathway-level transcriptional changes. Nominal or exploratory pathway patterns may still be useful for hypothesis generation, but they should not be reported as definitive pathway findings.
```

## Recommended Final Reporting Language

Use cautious language:

- "Topiramate treatment was associated with..."
- "The pathway-level results suggest..."
- "The strongest enrichment signal involved..."
- "The leading-edge genes indicate that the signal was driven by..."
- "These findings are exploratory unless they pass FDR < 0.05 or are replicated."

Avoid overclaiming:

- Do not write: "Topiramate caused activation of pathway X."
- Do not write: "Pathway X proves the mechanism of Topiramate."
- Do not write: "Genes in this pathway are biomarkers" unless a separate biomarker validation analysis was performed.

## Suggested Tables for the Report

### Table 1. Top GSEA Pathways

| Timepoint | Collection | Pathway | NES | FDR | Direction | Interpretation |
|---|---|---|---:|---:|---|---|
| Week 8 |  |  |  |  |  |  |
| Week 12 |  |  |  |  |  |  |

### Table 2. ORA-Supported Pathway Themes

| Timepoint | Threshold | Collection | Pathway | FDR | Overlap genes | Theme |
|---|---|---|---|---:|---|---|
| Week 8 |  |  |  |  |  |  |
| Week 12 |  |  |  |  |  |  |

### Table 3. Integrated Biological Themes

| Biological theme | Supporting pathways | Timepoint(s) | Direction | Driver genes | Strength of evidence |
|---|---|---|---|---|---|
|  |  |  |  |  |  |

## Short Conclusion Template

```text
The pathway analysis of GSE107015 suggests that Topiramate treatment is associated with coordinated whole-blood transcriptional changes related to [MAIN BIOLOGICAL THEMES]. The strongest evidence was observed at [WEEK 8 / WEEK 12 / BOTH TIMEPOINTS], where [PATHWAYS] showed [positive/negative] enrichment. ORA results [supported/did not support] these themes through overlap with [GENES/PATHWAYS]. These results should be interpreted as pathway-level treatment-associated signatures and should be considered exploratory unless supported by strict FDR control, consistency across timepoints, and independent validation.
```

## Interpretation Checklist

Before writing the final discussion, confirm that:

- [ ] GSEA results were interpreted using `padj`, not only raw `pval`.
- [ ] The direction of enrichment was based on the sign of `NES`.
- [ ] Similar pathways were grouped into themes.
- [ ] Leading-edge genes were reviewed for the main findings.
- [ ] ORA was used as supporting evidence, not as the only evidence.
- [ ] Week 8 and week 12 were compared separately.
- [ ] Conclusions were phrased as associations, not causal proof.
- [ ] Whole-blood tissue context was acknowledged.

## References

- GSE107015 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107015
- GSEA-MSigDB documentation: https://docs.gsea-msigdb.org/
- MSigDB human collections: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- MSigDB Hallmark collection details: https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp
- Reactome pathway database: https://reactome.org/
