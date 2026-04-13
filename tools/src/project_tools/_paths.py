from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Paths:
    root: Path
    data_dir: Path
    artifacts_dir: Path
    r_dir: Path
    raw_dir: Path
    r_library_dir: Path
    series_matrix_path: Path
    master_sample_sheet_path: Path
    expression_matrix_csv: Path
    preprocessing_qc_csv: Path

    @classmethod
    def from_pipeline_file(cls, path, level):
        root = Path(path).resolve().parents[level]
        data_dir = root / "data"
        artifacts_dir = root / "artifacts"
        r_dir = root / "R"
        return cls(
            root=root,
            data_dir=data_dir,
            artifacts_dir=artifacts_dir,
            r_dir=r_dir,
            raw_dir=data_dir / "GSE107015_RAW",
            r_library_dir=root / ".r_library",
            series_matrix_path=data_dir / "GSE107015_series_matrix.txt",
            master_sample_sheet_path=artifacts_dir / "master_sample_sheet.csv",
            expression_matrix_csv=artifacts_dir / "expression_matrix_rma.csv",
            preprocessing_qc_csv=artifacts_dir / "preprocessing_qc.csv",
        )


INTERNAL_CONFIGURATION = Paths.from_pipeline_file(__file__, 3)
