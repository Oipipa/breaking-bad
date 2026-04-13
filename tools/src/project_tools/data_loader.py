import gzip
import shutil
import tarfile
from urllib.request import urlretrieve
from project_tools._paths import INTERNAL_CONFIGURATION

SERIES_MATRIX_GZ_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/matrix/GSE107015_series_matrix.txt.gz"
RAW_TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/suppl/GSE107015_RAW.tar"


def _download_series_matrix_archive(path_configuration):
    series_matrix_archive = path_configuration.data_dir / "GSE107015_series_matrix.txt.gz"
    if not series_matrix_archive.exists():
        series_matrix_archive.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {series_matrix_archive.name}")
        urlretrieve(SERIES_MATRIX_GZ_URL, series_matrix_archive)
    return series_matrix_archive


def _extract_series_matrix(path_configuration):
    series_matrix_archive = _download_series_matrix_archive(path_configuration)
    series_matrix_path = path_configuration.data_dir / "GSE107015_series_matrix.txt"
    if not series_matrix_path.exists():
        print(f"Extracting {series_matrix_path.name}")
        with gzip.open(series_matrix_archive, "rb") as src, series_matrix_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return series_matrix_path


def _download_raw_archive(path_configuration):
    raw_archive = path_configuration.data_dir / "GSE107015_RAW.tar"
    if not raw_archive.exists():
        raw_archive.parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {raw_archive.name}")
        urlretrieve(RAW_TAR_URL, raw_archive)
    return raw_archive


def _extract_raw_data(path_configuration):
    raw_dir = path_configuration.raw_dir
    raw_archive = _download_raw_archive(path_configuration)
    if not any(raw_dir.glob("*.CEL.gz")):
        raw_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting {raw_archive.name}")
        with tarfile.open(raw_archive) as tar:
            tar.extractall(raw_dir)
    return raw_dir


def _remove_archives(path_configuration):
    (path_configuration.data_dir / "GSE107015_series_matrix.txt.gz").unlink(missing_ok=True)
    (path_configuration.data_dir / "GSE107015_RAW.tar").unlink(missing_ok=True)


def download_gse107015_data(path_configuration):
    matrix_txt = path_configuration.data_dir / "GSE107015_series_matrix.txt"
    raw_dir = path_configuration.raw_dir

    if matrix_txt.exists():
        print(f"Using existing {matrix_txt.name}")
    else:
        _extract_series_matrix(path_configuration)

    if any(raw_dir.glob("*.CEL.gz")):
        print(f"Using existing {raw_dir.name}/")
    else:
        _extract_raw_data(path_configuration)

    _remove_archives(path_configuration)

if __name__ == "__main__":
    download_gse107015_data(INTERNAL_CONFIGURATION)
