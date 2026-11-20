"""Small, explicit JSON/CSV contracts, with atomic output writes."""
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path


def read_config(path):
    path = Path(path).resolve()
    with path.open(encoding="utf-8-sig") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError("Configuration must be a mapping: {}".format(path))
    return value


def resolve_path(base, value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (Path(base) / path).resolve()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
        os.replace(tmp, str(path))
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def read_csv(path, required=()):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError("{} requires columns {}".format(path, ", ".join(required)))
        rows = []
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed CSV row {} in {}".format(reader.line_num, path))
            if any(v for v in row.values()):
                rows.append({k: v.strip() for k, v in row.items()})
        return rows


def write_csv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(value, name):
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError("{} must be a finite number, got {!r}".format(name, value)) from exc
    if not math.isfinite(result):
        raise ValueError("{} must be finite".format(name))
    return result


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
