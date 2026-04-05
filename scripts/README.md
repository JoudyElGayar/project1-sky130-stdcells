# Scripts

Three Python scripts automate the full characterization pipeline.

---

## characterize.py

Generates SPICE netlists, runs ngspice in batch mode, and extracts NLDM timing tables.

### Usage

```bash
cd project1/scripts

python3 characterize.py \
  --sky130 /path/to/sky130A \
  [--workdir ../netlists] \
  [--cells invx1 invx2 ...]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--sky130` | Yes | Path to sky130A directory (must contain `libs.ref/sky130_fd_pr/spice/`) |
| `--workdir` | No | Directory containing `stdcells.spice` and `sky130_params.spice` (default: `../netlists`) |
| `--cells` | No | Space-separated list of cells to run (default: all 13) |

### Output

- `results/<cell_name>/` — one `.sp` netlist and `.sp.log` per simulation point
- `results/<cell_name>.json` — 4 NLDM tables (7×7 arrays) for that cell
- `results/nldm_tables.json` — aggregated results for the batch that was run

### Dependencies

```bash
pip install numpy
# ngspice must be on your PATH
```

---

## merge_results.py

Merges per-cell JSON files produced by separate `characterize.py` runs into a single `nldm_tables.json`. Needed when cells were characterized in separate batches (since each `characterize.py` run overwrites `nldm_tables.json` with only the cells it processed).

### Usage

```bash
cd project1
python3 scripts/merge_results.py
```

Looks for `results/invx1.json`, `results/invx2.json`, etc. and combines them. Prints a warning for any missing cells.

---

## generate_report.py

Reads `results/nldm_tables.json` and produces plots and a markdown report.

### Usage

```bash
cd project1
python3 scripts/generate_report.py
```

### Output

| File | Description |
|------|-------------|
| `results/inv_delay_vs_load.pdf` | Cell_rise and cell_fall vs load for invx1–invx8 at tin=0.1225 ns |
| `results/rc_vs_spice.pdf` | RC model vs SPICE for invx1 over all 7 load values |
| `results/nldm_report.md` | Full NLDM tables for all 13 cells in markdown format |

Also prints to console:
- RC model parameters (Ron, Cg)
- All 13 × 4 NLDM tables formatted as grids
- RC vs SPICE % error at the midpoint condition (tin=0.1225 ns, Cload=0.0094 pF)

### Dependencies

```bash
pip install numpy matplotlib
```
