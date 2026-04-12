from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PipelinePaths:
    root: Path
    data_dir: Path
    artifacts_dir: Path
    r_dir: Path
    master_sheet_path: Path
    raw_dir: Path

    @classmethod
    def from_pipeline_file(cls, path: Path) -> "PipelinePaths":
        root = Path(path).resolve().parents[1]
        data_dir = root / "data"
        artifacts_dir = root / "artifacts"
        r_dir = root / "R"
        return cls(
            root=root,
            data_dir=data_dir,
            artifacts_dir=artifacts_dir,
            r_dir=r_dir,
            master_sheet_path=artifacts_dir / "master_sample_sheet.csv",
            raw_dir=data_dir / "GSE107015_RAW"
        )


PATHS = PipelinePaths.from_pipeline_file(__file__)
