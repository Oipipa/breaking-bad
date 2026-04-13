from project_tools._paths import Paths
from project_tools.data_loader import download_gse107015_data
from project_tools.sample_sheet_generator import build_master_sample_sheet 
from project_tools.expression_preprocessing import run_rma_preprocessing
from project_tools.preprocessing_results import quality_control_summary

PATH_CONFIGURATION = Paths.from_pipeline_file(__file__, 1)

def execute_pipeline(): 
    download_gse107015_data(PATH_CONFIGURATION)
    build_master_sample_sheet(PATH_CONFIGURATION)
    run_rma_preprocessing(PATH_CONFIGURATION)
    quality_control_summary(PATH_CONFIGURATION)

if __name__ == "__main__": 
    execute_pipeline()