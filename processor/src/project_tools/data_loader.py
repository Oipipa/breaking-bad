import gzip
import shutil
import tarfile
from pathlib import Path
from urllib.request import urlretrieve

SERIES_MATRIX_GZ_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/matrix/GSE107015_series_matrix.txt.gz"
RAW_TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/suppl/GSE107015_RAW.tar"


def _download(url: str, path: Path) -> Path:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {path.name}")
        urlretrieve(url, path)
    return path


def _gunzip(source: Path, target: Path) -> Path:
    if not target.exists():
        print(f"Extracting {target.name}")
        with gzip.open(source, "rb") as src, target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return target


def _untar(source: Path, target_dir: Path) -> Path:
    if not any(target_dir.glob("*.CEL.gz")):
        target_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting {source.name}")
        with tarfile.open(source) as tar:
            tar.extractall(target_dir)
    return target_dir


def download_gse107015_data(data_dir, *, keep_archives: bool = False) -> dict[str, Path]:
    matrix_txt = data_dir / "GSE107015_series_matrix.txt"
    raw_dir = data_dir / "GSE107015_RAW"
    matrix_gz = data_dir / "GSE107015_series_matrix.txt.gz"
    raw_tar = data_dir / "GSE107015_RAW.tar"

    if matrix_txt.exists():
        print(f"Using existing {matrix_txt.name}")
    else:
        _gunzip(_download(SERIES_MATRIX_GZ_URL, matrix_gz), matrix_txt)

    if any(raw_dir.glob("*.CEL.gz")):
        print(f"Using existing {raw_dir.name}/")
    else:
        _untar(_download(RAW_TAR_URL, raw_tar), raw_dir)

    if not keep_archives:
        matrix_gz.unlink(missing_ok=True)
        raw_tar.unlink(missing_ok=True)

    return {"data_dir": data_dir, "series_matrix": matrix_txt, "raw_dir": raw_dir}
