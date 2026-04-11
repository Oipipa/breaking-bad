import csv
import gzip
import urllib.request
from pathlib import Path
import re

SERIES_MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107015/matrix/GSE107015_series_matrix.txt.gz"
FILENAME_PATTERN = re.compile(
    r"^(?P<geo_sample_id>GSM\d+)_"
    r"(?P<batch_code>EA05064)_"
    r"(?P<array_scan_id>\d+)"
    r"(?:_(?P<filename_flag>[A-Z]\d))?_"
    r"(?P<platform_code>H133_Plus)_"
    r"(?P<subject_id>\d+)_"
    r"(?P<timepoint_code>[234])\.CEL\.gz$"
)
TIMEPOINT_MAP = {"Baseline": "baseline", "Week 8": "week8", "Week 12": "week12"}
CHARACTERISTIC_COLUMNS = {
    "individual id": "geo_subject_id",
    "gender": "gender",
    "age": "age",
    "race": "race",
    "ethnicity": "ethnicity",
    "treatment": "treatment_arm",
    "time point": "geo_timepoint",
    "tissue": "tissue"
}

def _read_series_matrix():
    request = urllib.request.Request(SERIES_MATRIX_URL, headers={"User-Agent": "python"})
    payload = urllib.request.urlopen(request).read()
    text = gzip.decompress(payload).decode("utf-8", errors="replace")
    sample_rows = {}
    characteristic_rows = []
    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue
        row = next(csv.reader([line], delimiter="\t"))
        key, values = row[0], row[1:]
        if key == "!Sample_characteristics_ch1":
            characteristic_rows.append(values)
        elif key not in sample_rows:
            sample_rows[key] = values
    geo_ids = sample_rows["!Sample_geo_accession"]
    metadata = {geo_id: {} for geo_id in geo_ids}
    direct_columns = {"!Sample_title": "sample_title", "!Sample_source_name_ch1": "source_name", "!Sample_platform_id": "platform_id"}
    for key, column in direct_columns.items():
        for geo_id, value in zip(geo_ids, sample_rows[key]):
            metadata[geo_id][column] = value
    for values in characteristic_rows:
        label = values[0].split(":", 1)[0].strip().lower()
        column = CHARACTERISTIC_COLUMNS.get(label)
        if column is None:
            continue
        for geo_id, value in zip(geo_ids, values):
            metadata[geo_id][column] = value.split(":", 1)[1].strip()
    return metadata


def build_master_sample_sheet(data_dir, output_path):
    geo_metadata = _read_series_matrix()
    rows = []
    for raw_file in sorted(data_dir.glob("*.CEL.gz")):
        match = FILENAME_PATTERN.match(raw_file.name)
        row = match.groupdict()
        geo_row = geo_metadata[row["geo_sample_id"]]
        row["series_id"] = "GSE107015"
        row["sample_title"] = geo_row["sample_title"]
        row["source_name"] = geo_row["source_name"]
        row["subject_id"] = geo_row["geo_subject_id"]
        row["treatment_arm"] = geo_row["treatment_arm"]
        row["timepoint"] = TIMEPOINT_MAP[geo_row["geo_timepoint"]]
        row["gender"] = geo_row["gender"]
        row["age"] = geo_row["age"]
        row["race"] = geo_row["race"]
        row["ethnicity"] = geo_row["ethnicity"]
        row["tissue"] = geo_row["tissue"]
        row["platform_id"] = geo_row["platform_id"]
        row["raw_filename"] = raw_file.name
        row["cel_filename"] = raw_file.name.removesuffix(".gz")
        row.pop("filename_flag", None)
        rows.append(row)
    rows.sort(key=lambda row: (int(row["subject_id"]), int(row["timepoint_code"]), row["geo_sample_id"]))
    fieldnames = [
        "series_id",
        "geo_sample_id",
        "sample_title",
        "source_name",
        "subject_id",
        "treatment_arm",
        "timepoint",
        "timepoint_code",
        "raw_filename",
        "cel_filename",
        "gender",
        "age",
        "race",
        "ethnicity",
        "tissue",
        "batch_code",
        "array_scan_id",
        "platform_id",
        "platform_code"
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)