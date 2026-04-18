This project is the implementation part of the Big Data for Bioanalytics and Medicine (S2_BDBM6ILV) course taught in the summer semester of 2026 at IMC FH Krems.

The data source for this project is [GSE107015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107015) and the series matrix can be downloaded [here](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/matrix/GSE107015_series_matrix.txt.gz). This project will be focused on reproducing the results in [this article](https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2017.00271/full) titled "Identification of Novel Signal Transduction, Immune Function, and Oxidative Stress Genes and Pathways by Topiramate for Treatment of Methamphetamine Dependence Based on Secondary Outcomes". 

## Project Structure

- `tools/`: the Python project. This contains the installable package for any reusable utilities for the pipeline.
- `pipeline/`: the Python files here use the logic defined in `tools/` to illustrate the workflow in a more condensed view. The files contained here do not implement the actual analysis logic and mainly define parameters and call pre-built functions.
- `R/`: the R part of the project. This contains standalone R scripts for steps that are easier or more appropriate to run in R.
- `data/` (not committed to repo): raw input data. Downloaded source files, raw archives, and any large local-only inputs should be here.
- `artifacts/`: derived outputs. Preprocessed tables, metadata products, cached intermediate results, and other generated artifacts should be here.

## Data Conventions

The intended convention is to keep raw data and generated data separate.

- We keep raw or externally downloaded inputs in `data/` at the project's root (this is gitignored). **The pipeline downloads and extract the data automatically, so manual downloads are not required**. 
- We keep preprocessed outputs and other generated artifacts in `artifacts/` (this is gitignored because artifacts like the expression matrix are simply too large to commit).
- We treat `data/` as the source location for inputs and `artifacts/` as the destination for anything produced by the pipeline.

## Python Environments

This project was scaffolded with `uv` in mind, so it is the recommended tool for managing your virtual environments, though you can use anything else that works (eg. conda).

So, the intended workflow is: 

```bash
uv venv 
source .venv/bin/activate
cd tools
uv pip install -e .
```

Once `tools/` is installed, notebooks (or external scripts) can import from the package directly, assuming they use the same environment.

Example:

```python
from project_tools.hello_world import hello_world

hello_world("Anubhav")
```

## Running the pipeline 

The pipeline itself is a python package at the root of the directory, so it can be executed as: 

```bash
python -m pipeline
```

## limma Differential Expression (Per Week)

The pipeline includes a per-week limma workflow on top of the RMA output.

The DE model is run on baseline-normalized expression deltas:

- week 8 delta = expression(week8) - expression(baseline)
- week 12 delta = expression(week12) - expression(baseline)

### Analysis mode

The workflow runs in public no-label mode only.

At each timepoint, it compares Topiramate vs Placebo on baseline-normalized deltas:

- week 8: `Topiramate - Placebo` on `(week8 - baseline)`
- week 12: `Topiramate - Placebo` on `(week12 - baseline)`

### limma outputs

Outputs are written to `artifacts/limma/`.

Public no-label mode outputs:

- `de_public_week8_topiramate_vs_placebo_limma.csv`
- `de_public_week12_topiramate_vs_placebo_limma.csv`
- `overlap_public_week8_week12.csv`
- `overlap_public_week8_week12_pathway_input.csv`
- `limma_comparison_summary.csv`
- `limma_metadata_summary.csv`

Pathway-ready files (generated per DE table):

- `*_pathway_input_all.csv`
- `*_pathway_input_p05.csv`
- `*_pathway_input_fdr015.csv`
- `*_ranked.rnk`

## How to Read the Results

- `limma_comparison_summary.csv` gives gene-level signal counts per contrast.
