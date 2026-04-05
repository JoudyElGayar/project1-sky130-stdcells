
"""
merge_results.py
Reads all per-cell JSON files from results/ and merges them into nldm_tables.json
Run from the project1/ directory: python3 scripts/merge_results.py
"""
import json
from pathlib import Path

results_dir = Path("results")
cells = [
    "invx1","invx2","invx4","invx8",
    "nand2x1","nand2x2","nand2x4",
    "nor2x1","nor2x2","nor2x4",
    "maj3x1","maj3x2","maj3x4",
]

all_results = {}
missing = []

for cell in cells:
    p = results_dir / f"{cell}.json"
    if p.exists():
        with open(p) as f:
            all_results[cell] = json.load(f)
        print(f"  loaded {cell}.json")
    else:
        missing.append(cell)
        print(f"  MISSING: {cell}.json")

if missing:
    print(f"\nWARNING: {len(missing)} cell(s) missing: {missing}")
    print("Re-run characterize.py for those cells first.")
else:
    print(f"\nAll 13 cells loaded.")

out = results_dir / "nldm_tables.json"
with open(out, "w") as f:
    json.dump(all_results, f, indent=2)
print(f"Saved {len(all_results)} cells -> {out}")
