"""Build a source archive from an explicit set of project files."""
import hashlib
import os
import zipfile
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]
VERSION = "0.2.0"
DIRECTORIES = ("src", "docs", "examples", "tests", "tools")
FILES = ("README.md", "NOTICE.md", "setup.py", "setup.cfg", "pyproject.toml", "MANIFEST.in", ".gitignore")


def main():
    destination = PACKAGE / "dist" / ("CO2RR_Network_" + VERSION + "_source.zip")
    destination.parent.mkdir(parents=True, exist_ok=True)
    root = Path("CO2RR_Network_" + VERSION)
    with zipfile.ZipFile(str(destination), "w", zipfile.ZIP_DEFLATED) as archive:
        for filename in FILES:
            archive.write(str(PACKAGE / filename), str(root / filename))
        for name in DIRECTORIES:
            for directory, children, filenames in os.walk(str(PACKAGE / name)):
                children[:] = sorted(c for c in children if c not in ("__pycache__", "fine-tuning") and not c.endswith(".egg-info"))
                for filename in sorted(filenames):
                    path = Path(directory) / filename
                    if path.suffix not in (".pyc", ".pyo"):
                        archive.write(str(path), str(root / path.relative_to(PACKAGE)))
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    destination.with_suffix(".zip.sha256").write_text(digest + "  " + destination.name + "\n", encoding="utf-8")
    print(str(destination))


if __name__ == "__main__":
    main()
