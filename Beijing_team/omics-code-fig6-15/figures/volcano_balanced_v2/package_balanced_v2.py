"""Package the two revised volcano PDFs with the other thirteen vector panels."""

from __future__ import annotations

from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ORIGINAL_PANELS = PARENT / "panels"
REPLACED = {
    "Fig6_J_volcano_48pt_Arial.pdf": HERE / "Fig6_J_volcano_Arial48_balanced_v2.pdf",
    "FigS2_A_volcano_48pt_Arial.pdf": HERE / "FigS2_A_volcano_Arial48_balanced_v2.pdf",
}

all_vector = sorted(
    p for p in ORIGINAL_PANELS.glob("*.pdf")
    if "_original_art" not in p.name
)
assert len(all_vector) == 15, [p.name for p in all_vector]
assert all(p.is_file() for p in REPLACED.values())

target = HERE / "Fig6_JM_FigS2_15_vector_PDFs_Arial48_balanced_v2.zip"
with ZipFile(target, "w", ZIP_DEFLATED, compresslevel=8) as archive:
    for source in all_vector:
        actual = REPLACED.get(source.name, source)
        archive.write(actual, f"panels/{source.name}")
    archive.write(HERE / "README.md", "README.md")
    archive.write(HERE / "draw_volcano_balanced_v2.py", "draw_volcano_balanced_v2.py")

with ZipFile(target) as archive:
    assert archive.testzip() is None
    assert len([n for n in archive.namelist() if n.endswith(".pdf")]) == 15
print(target)
