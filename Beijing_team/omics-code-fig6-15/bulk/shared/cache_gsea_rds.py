"""Convert the two saved R GSEA objects into a Python pickle cache."""

from pathlib import Path
import pickle
import rdata


ROOT = Path(r"E:\PCP\bulk\20260803\turn_20260802065502_636832")
OUT = ROOT / "_spreadsheet_work" / "gsea_selection" / "gsea_rds_cache.pkl"
files = {
    "PDC_vs_PDO": ROOT / "09_GSEA_visualization" / "gsea_objects_PDO_vs_PDC.RDS",
    "Ulixertinib_vs_Ctrl": ROOT / "09_GSEA_visualization" / "gsea_objects_Ulixertinib_vs_Ctrl.RDS",
}

cache = {}
for comparison, path in files.items():
    print(f"Parsing {comparison}: {path.name}", flush=True)
    cache[comparison] = rdata.conversion.convert(rdata.parser.parse_file(path))
    for key, value in cache[comparison].items():
        if hasattr(value, "__dict__"):
            fields = {k: type(v).__name__ for k, v in value.__dict__.items()}
            print(comparison, key, fields, flush=True)
        else:
            print(comparison, key, type(value).__name__, getattr(value, "shape", None), flush=True)

with OUT.open("wb") as handle:
    pickle.dump(cache, handle, protocol=pickle.HIGHEST_PROTOCOL)
print(f"Saved {OUT} ({OUT.stat().st_size} bytes)")
