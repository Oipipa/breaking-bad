## Premblers

These are the two helper modules used by `pipeline/preamble.py`.

- `src/project_tools/data_loader.py`: downloads the GSE107015 series matrix and raw CEL archive from GEO, extracts them into `data/`, reuses files if they already exist, and removes the downloaded archives by default.
- `src/project_tools/sample_sheet_generator.py`: reads the series matrix and the raw CEL filenames, pulls out the sample metadata needed for the project, and writes a master sample sheet CSV for the downstream analysis steps.

### Raw Data and Sample Sheet

This preamble uses the public GEO release for GSE107015. That matches the paper at a high level. 

- `GSE107015_series_matrix.txt`: the GEO series matrix. Here we only use the `!Sample_*` metadata rows.
- `GSE107015_RAW/*.CEL.gz`: the raw Affymetrix Human Genome U-133 Plus 2.0 files. Filenames such as `GSM2859594_EA05064_31476_H133_Plus_10265_3.CEL.gz` carry the GEO sample id, batch code, scan id, platform code, subject id, and a timepoint code.

That lines up with the paper's trial structure: whole-blood samples from placebo and topiramate arms at baseline, week 8, and week 12. The main difference is that this step follows the public GEO submission, so it works from the released array files and sample annotations, not the full internal set of blood draws or the responder/non-responder labels used later in the paper's secondary analysis.

The master sample sheet is just a join between the CEL filenames and the series-matrix metadata. It keeps one row per array with subject id, treatment arm, timepoint, demographics, platform and batch fields, and the original filenames so downstream steps can work with the same trial layout as the study.