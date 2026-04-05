"""
SKY130 Standard Cell NLDM Characterization Script
Generates SPICE netlists, runs ngspice, and extracts timing tables.

Output:
    results/  -> per-cell JSON data
    nldm_tables.json -> aggregated NLDM tables for all 13 cells
"""

import os
import re
import json
import subprocess
import argparse
import itertools
import numpy as np
from pathlib import Path

# ─────────────────────────────────────────────
# NLDM sweep vectors (from project spec)
# ─────────────────────────────────────────────
INPUT_TRANSITIONS_NS = [0.0100, 0.0231, 0.0531, 0.1225, 0.2823, 0.6507, 1.5000]
LOAD_CAPS_PF         = [0.0005, 0.0013, 0.0035, 0.0094, 0.0249, 0.0662, 0.1758]

VDD      = 1.8          # Supply voltage
VTH      = VDD * 0.5   # Propagation delay threshold  (0.9 V)
V20      = VDD * 0.2   # Slew low  threshold          (0.36 V)
V80      = VDD * 0.8   # Slew high threshold           (1.44 V)
TEMP     = 25
CORNER   = "tt"

# ─────────────────────────────────────────────
# Cell definitions
#   name : (subckt_name, input_pin, output_pin, static_inputs)
#   static_inputs: dict of extra pin=voltage needed to activate the path
# ─────────────────────────────────────────────
# For multi-input cells we pick the worst-case / standard input to toggle.
# NAND2: toggle A, B=1 (so pull-down stack can conduct)
# NOR2 : toggle A, B=0 (so pull-up stack can conduct)
# MAJ3 : toggle A, B=1 C=1 (majority changes with A)
CELLS = {
    # name          subckt     toggle  out   static pins (pin: V)
    "invx1":  ("invx1",  "A", "Z", {}),
    "invx2":  ("invx2",  "A", "Z", {}),
    "invx4":  ("invx4",  "A", "Z", {}),
    "invx8":  ("invx8",  "A", "Z", {}),
    "nand2x1":("nand2x1","A", "Z", {"B": VDD}),
    "nand2x2":("nand2x2","A", "Z", {"B": VDD}),
    "nand2x4":("nand2x4","A", "Z", {"B": VDD}),
    "nor2x1": ("nor2x1", "A", "Z", {"B": 0.0}),
    "nor2x2": ("nor2x2", "A", "Z", {"B": 0.0}),
    "nor2x4": ("nor2x4", "A", "Z", {"B": 0.0}),
    "maj3x1": ("maj3x1", "A", "Z", {"B": VDD, "C": 0.0}),
    "maj3x2": ("maj3x2", "A", "Z", {"B": VDD, "C": 0.0}),
    "maj3x4": ("maj3x4", "A", "Z", {"B": VDD, "C": 0.0}),
}

# Port order for each subckt (needed to build Xcell line)
SUBCKT_PORTS = {
    "invx1":   ["A", "Z", "VDD", "GND"],
    "invx2":   ["A", "Z", "VDD", "GND"],
    "invx4":   ["A", "Z", "VDD", "GND"],
    "invx8":   ["A", "Z", "VDD", "GND"],
    "nand2x1": ["A", "B", "Z", "VDD", "GND"],
    "nand2x2": ["A", "B", "Z", "VDD", "GND"],
    "nand2x4": ["A", "B", "Z", "VDD", "GND"],
    "nor2x1":  ["A", "B", "Z", "VDD", "GND"],
    "nor2x2":  ["A", "B", "Z", "VDD", "GND"],
    "nor2x4":  ["A", "B", "Z", "VDD", "GND"],
    "maj3x1":  ["A", "B", "C", "Z", "VDD", "GND"],
    "maj3x2":  ["A", "B", "C", "Z", "VDD", "GND"],
    "maj3x4":  ["A", "B", "C", "Z", "VDD", "GND"],
}


def build_netlist(cell_name, subckt, toggle_pin, out_pin, static_pins,
                  tin_ns, cload_pf, sky130_dir, work_dir):
    """
    Generate a transient SPICE netlist for one (cell, tin, cload) point.
    Returns the netlist string and a dict of measure names.
    """
    ports   = SUBCKT_PORTS[subckt]
    tin_s   = tin_ns  * 1e-9
    cload_f = cload_pf * 1e-12

    # Simulation time: enough for two full transitions + settling
    Ron_est     = 3000
    RC_ns       = Ron_est * (cload_pf * 1e-12) * 1e9
    sim_time_ns = max(20.0, 10 * tin_ns + 80 * RC_ns)
    sim_time_ns = min(sim_time_ns, 800.0)
    tstep_ns    = sim_time_ns / 20000

    td  = 3 * tin_s
    pw  = max(10 * tin_s, 20 * RC_ns * 1e-9)
    per = 2 * (td + pw + 5 * RC_ns * 1e-9)
    if per * 1e9 > sim_time_ns * 0.9:
        sim_time_ns = min(per * 1e9 / 0.9, 800.0)

    nfet_model = os.path.join(
        sky130_dir,
        "libs.ref/sky130_fd_pr/spice/sky130_fd_pr__nfet_01v8__tt.pm3.spice")
    pfet_model = os.path.join(
        sky130_dir,
        "libs.ref/sky130_fd_pr/spice/sky130_fd_pr__pfet_01v8__tt.pm3.spice")

    lines = []
    lines.append(f"* Characterization: {cell_name}  tin={tin_ns}ns  cload={cload_pf}pF")
    lines.append(f".temp {TEMP}")
    lines.append(f'.include "{os.path.join(work_dir, "sky130_params.spice")}"')
    lines.append(f'.include "{nfet_model}"')
    lines.append(f'.include "{pfet_model}"')
    lines.append(f'.include "{os.path.join(work_dir, "stdcells.spice")}"')
    lines.append("")

    # Supply
    lines.append(f"VVDD VDD GND {VDD}")

    static_node = {}
    for pin, vol in static_pins.items():
        node = f"static_{pin}"
        lines.append(f"V{pin} {node} GND {vol}")
        static_node[pin] = node

    # Toggling pulse source on toggle_pin
    lines.append(
        f"VIN {toggle_pin} GND "
        f"PULSE(0 {VDD} {td:.6e} {tin_s:.6e} {tin_s:.6e} {pw:.6e} {per:.6e})"
    )

    # Load cap
    lines.append(f"CL {out_pin} GND {cload_f:.6e}")

    # Cell instantiation
    node_map = {}
    for p in ports:
        if p == "VDD":
            node_map[p] = "VDD"
        elif p == "GND":
            node_map[p] = "GND"
        elif p == out_pin:
            node_map[p] = out_pin
        elif p == toggle_pin:
            node_map[p] = toggle_pin
        elif p in static_node:
            node_map[p] = static_node[p]
        else:
            node_map[p] = p
    node_list = " ".join(node_map[p] for p in ports)
    lines.append(f"XCELL {node_list} {subckt}")
    lines.append("")

    # Transient
    lines.append(f".tran {tstep_ns:.6e}n {sim_time_ns:.6e}n")
    lines.append("")

    # ── Measurements ──────────────────────────────────────────────
    # We use FALL=1/RISE=1 referencing the FIRST crossing after t=0
    # cell_rise : input falls (INV: A falls -> Z rises), so TRIG A VAL VTH FALL=1
    # cell_fall : input rises (INV: A rises -> Z falls), so TRIG A VAL VTH RISE=1
    # For cells where output is non-inverting (MAJ3 when B=C=1: A rises -> Z rises)
    # we swap. Detect inversion by cell type.

    vthr = round(VTH, 4)   # 0.9
    v20r = round(V20, 4)   # 0.36
    v80r = round(V80, 4)   # 1.44

    inverting = subckt.startswith("inv") or subckt.startswith("nand") or subckt.startswith("nor") or subckt.startswith("maj")

    if inverting:
        lines.append(
            f".measure TRAN cell_rise "
            f"TRIG v({toggle_pin}) VAL={vthr} FALL=1 "
            f"TARG v({out_pin}) VAL={vthr} RISE=1"
        )
        lines.append(
            f".measure TRAN cell_fall "
            f"TRIG v({toggle_pin}) VAL={vthr} RISE=1 "
            f"TARG v({out_pin}) VAL={vthr} FALL=1"
        )
    else:
        lines.append(
            f".measure TRAN cell_rise "
            f"TRIG v({toggle_pin}) VAL={vthr} RISE=1 "
            f"TARG v({out_pin}) VAL={vthr} RISE=1"
        )
        lines.append(
            f".measure TRAN cell_fall "
            f"TRIG v({toggle_pin}) VAL={vthr} FALL=1 "
            f"TARG v({out_pin}) VAL={vthr} FALL=1"
        )
    lines.append(
        f".measure TRAN rise_transition "
        f"TRIG v({out_pin}) VAL={v20r} RISE=1 "
        f"TARG v({out_pin}) VAL={v80r} RISE=1"
    )
    lines.append(
        f".measure TRAN fall_transition "
        f"TRIG v({out_pin}) VAL={v80r} FALL=1 "
        f"TARG v({out_pin}) VAL={v20r} FALL=1"
    )
    lines.append("")
    lines.append(".end")
    return "\n".join(lines)


def run_ngspice(netlist_path):
    """Run ngspice in batch mode, return stdout."""
    result = subprocess.run(
        ["ngspice", "-b", "-o", str(netlist_path) + ".log", str(netlist_path)],
        capture_output=True, text=True, timeout=120
    )
    return result.stdout + result.stderr


def parse_measure(output, measure_name):
    """
    Extract a .measure result value from ngspice output.
    Returns value in seconds (as printed by ngspice), or None.
    """
    # ngspice prints: measure_name = value
    pattern = re.compile(
        r"^\s*" + re.escape(measure_name) + r"\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
        re.IGNORECASE | re.MULTILINE
    )
    # Also try the log file format
    m = pattern.search(output)
    if m:
        val = float(m.group(1))
        return val if val > 0 else None
    return None


def parse_log_file(log_path, measure_names):
    """Parse the .log file ngspice writes."""
    try:
        with open(log_path, "r") as f:
            text = f.read()
        results = {}
        for name in measure_names:
            results[name] = parse_measure(text, name)
        return results
    except FileNotFoundError:
        return {n: None for n in measure_names}


def characterize_cell(cell_name, sky130_dir, work_dir, out_dir):
    """
    Run all 7×7 simulations for a cell.
    Returns dict with 4 tables (7×7 numpy arrays) in ns.
    """
    subckt, toggle_pin, out_pin, static_pins = CELLS[cell_name]

    tables = {
        "cell_rise":        np.full((7, 7), np.nan),
        "cell_fall":        np.full((7, 7), np.nan),
        "rise_transition":  np.full((7, 7), np.nan),
        "fall_transition":  np.full((7, 7), np.nan),
    }

    measure_names = ["cell_rise", "cell_fall", "rise_transition", "fall_transition"]

    total = len(INPUT_TRANSITIONS_NS) * len(LOAD_CAPS_PF)
    done  = 0

    for i, tin in enumerate(INPUT_TRANSITIONS_NS):
        for j, cload in enumerate(LOAD_CAPS_PF):
            done += 1
            netlist = build_netlist(
                cell_name, subckt, toggle_pin, out_pin, static_pins,
                tin, cload, sky130_dir, work_dir
            )
            sp_path = out_dir / f"{cell_name}_tin{i}_cload{j}.sp"
            with open(sp_path, "w") as f:
                f.write(netlist)

            print(f"  [{done:2d}/{total}] {cell_name}  tin={tin:.4f}ns  cload={cload:.4f}pF  ...",
                  end=" ", flush=True)

            try:
                raw_out = run_ngspice(sp_path)
            except subprocess.TimeoutExpired:
                print("TIMEOUT")
                continue
            except FileNotFoundError:
                print("ERROR: ngspice not found. Is it in your PATH?")
                raise

            log_path = str(sp_path) + ".log"
            vals = parse_log_file(log_path, measure_names)

            
            for name in measure_names:
                if vals[name] is None:
                    vals[name] = parse_measure(raw_out, name)

            # Convert seconds -> nanoseconds
            ok_count = 0
            for name in measure_names:
                v = vals[name]
                if v is not None:
                    row, col = i, j
                    tables[name][row, col] = v * 1e9   # s -> ns
                    ok_count += 1

            print(f"OK ({ok_count}/4 measures)")

    return tables


def fmt_table(arr, row_idx, col_idx):
    """Format a 2D numpy array as a nested list (JSON-serializable)."""
    return arr.tolist()


def main():
    parser = argparse.ArgumentParser(description="SKY130 NLDM Characterizer")
    parser.add_argument(
        "--sky130",
        required=True,
        help="Path to sky130A directory (contains libs.ref/sky130_fd_pr/spice/)")
    parser.add_argument(
        "--cells",
        nargs="+",
        default=list(CELLS.keys()),
        help="Subset of cells to characterize (default: all 13)")
    parser.add_argument(
        "--workdir",
        default=str(Path(__file__).resolve().parent.parent / "netlists"),
        help="Directory containing stdcells.spice and sky130_params.spice")
    args = parser.parse_args()

    sky130_dir = os.path.abspath(args.sky130)
    work_dir   = os.path.abspath(args.workdir)
    out_dir    = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True)

    # Verify required files exist
    for fname in ["stdcells.spice", "sky130_params.spice"]:
        p = os.path.join(work_dir, fname)
        if not os.path.exists(p):
            print(f"ERROR: {fname} not found in {work_dir}")
            return

    nfet = os.path.join(sky130_dir,
        "libs.ref/sky130_fd_pr/spice/sky130_fd_pr__nfet_01v8__tt.pm3.spice")
    pfet = os.path.join(sky130_dir,
        "libs.ref/sky130_fd_pr/spice/sky130_fd_pr__pfet_01v8__tt.pm3.spice")
    for p in [nfet, pfet]:
        if not os.path.exists(p):
            print(f"ERROR: Model file not found:\n  {p}")
            return

    all_results = {}

    for cell_name in args.cells:
        if cell_name not in CELLS:
            print(f"WARNING: Unknown cell '{cell_name}', skipping.")
            continue
        print(f"\n{'='*60}")
        print(f"Characterizing: {cell_name}")
        print(f"{'='*60}")

        cell_out_dir = out_dir / cell_name
        cell_out_dir.mkdir(exist_ok=True)

        tables = characterize_cell(cell_name, sky130_dir, work_dir, cell_out_dir)

        # Save per-cell JSON
        cell_data = {
            "cell": cell_name,
            "index1": INPUT_TRANSITIONS_NS,   # input_transition (ns)
            "index2": LOAD_CAPS_PF,            # capacitance (pF)
            "tables": {k: fmt_table(v, INPUT_TRANSITIONS_NS, LOAD_CAPS_PF)
                       for k, v in tables.items()}
        }
        with open(out_dir / f"{cell_name}.json", "w") as f:
            json.dump(cell_data, f, indent=2)

        all_results[cell_name] = cell_data
        print(f"  -> Saved {cell_name}.json")

    # Save aggregated results
    with open(out_dir / "nldm_tables.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n✓ All done. Results in: {out_dir}/nldm_tables.json")


if __name__ == "__main__":
    main()
