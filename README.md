# Project 1 — Standard Cell Library Design and Characterization (SKY130)

**Course:** VLSI Design  
**Process:** SkyWater SKY130 (180 nm)  
**Corner:** TT | **VDD:** 1.8 V | **Temperature:** 25°C

---

## Repository Structure

```
project1/
├── netlists/
│   ├── stdcells.spice        # SPICE subcircuit definitions for all 13 cells
│   ├── sky130_params.spice   # SKY130 TT-corner mismatch parameters
│   └── test_inv_fixed.spice  # Single-point test netlist for invx1 (verification)
│
├── scripts/
│   ├── characterize.py       # Main automation: generates netlists, runs ngspice, extracts tables
│   ├── generate_report.py    # Post-processing: plots, RC comparison, markdown tables
│   └── merge_results.py      # Utility: merges per-cell JSON files into nldm_tables.json
│
├── results/
│   ├── nldm_tables.json      # Aggregated NLDM tables for all 13 cells
│   ├── nldm_report.md        # Full NLDM table printout in markdown format
│   ├── inv_delay_vs_load.pdf # Plot: delay vs load for inverter family at tin=0.1225ns
│   ├── rc_vs_spice.pdf       # Plot: RC model vs SPICE comparison for invx1
│   ├── invx1.json            # Per-cell characterization data (one file per cell)
│   ├── invx2.json
│   └── ...                   # (invx4, invx8, nand2x1/2/4, nor2x1/2/4, maj3x1/2/4)
│
└── Project1_Report.pdf       # Final report (PDF)
```

---

## Cell Library

13 cells across 4 logic families:

| Family | Cells | Function |
|--------|-------|----------|
| Inverter | invx1, invx2, invx4, invx8 | Z = A' |
| NAND2 | nand2x1, nand2x2, nand2x4 | Z = (AB)' |
| NOR2 | nor2x1, nor2x2, nor2x4 | Z = (A+B)' |
| MAJ3 | maj3x1, maj3x2, maj3x4 | Z = AB + BC + AC |

Drive strengths scale transistor widths linearly (x2 = 2× widths, x4 = 4×, x8 = 8×).

---

## Reproducing the Results

### Requirements

- Python
- ngspice 
- numpy, matplotlib
- SKY130A PDK installed

### Step 1 — Run characterization

```bash
cd project1/scripts

python3 characterize.py \
  --sky130 /path/to/sky130A \
  --workdir ../netlists
```

This runs all 637 simulations (13 cells × 49 points each) and saves per-cell JSON files to `results/`.

To run a subset only:
```bash
python3 characterize.py --sky130 /path/to/sky130A --cells invx1 invx2
```

### Step 2 — Merge results (if run in batches)

```bash
cd project1
python3 scripts/merge_results.py
```

This combines all per-cell `.json` files into the aggregated `results/nldm_tables.json`.

### Step 3 — Generate report assets

```bash
cd project1
python3 scripts/generate_report.py
```

Outputs:
- `results/inv_delay_vs_load.pdf` — delay vs load plot for inverter family
- `results/rc_vs_spice.pdf` — RC model vs SPICE comparison
- `results/nldm_report.md` — full NLDM tables for all 13 cells

---

## NLDM Table Parameters

| Parameter | Values |
|-----------|--------|
| Input transitions (ns) | 0.0100, 0.0231, 0.0531, 0.1225, 0.2823, 0.6507, 1.5000 |
| Load capacitances (pF) | 0.0005, 0.0013, 0.0035, 0.0094, 0.0249, 0.0662, 0.1758 |
| Propagation delay threshold | 50% VDD = 0.9 V |
| Slew thresholds | 20%–80% VDD = 0.36 V – 1.44 V |
