from project_tools.expression_preprocessing import run_rma_preprocess
from __paths import PATHS

run_rma_preprocess(PATHS.raw_dir, PATHS.artifacts_dir, PATHS.r_dir)