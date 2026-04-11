from pathlib import Path
from project_tools.data_loader import download_gse107015_data
from project_tools.sample_sheet_generator import build_master_sample_sheet

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_CSV = ROOT / "artifacts" / "master_sample_sheet.csv"

download_gse107015_data(DATA_DIR)
build_master_sample_sheet(DATA_DIR, OUTPUT_CSV)