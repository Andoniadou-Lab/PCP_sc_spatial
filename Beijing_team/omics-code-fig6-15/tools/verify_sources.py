"""Verify collected source files against SOURCE_MANIFEST.csv (stdlib only)."""

import csv
import hashlib
from pathlib import Path

repo = Path(__file__).resolve().parents[1]
with (repo / "SOURCE_MANIFEST.csv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
bad = []
for row in rows:
    path = repo / row["repository_path"]
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        bad.append(row["repository_path"])
print(f"Source files verified: {len(rows) - len(bad)}/{len(rows)}")
for path in bad:
    print("Mismatch:", path)
raise SystemExit(bool(bad))
