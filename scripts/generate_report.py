#!/usr/bin/env python3
"""
generate_report.py
Reads nldm_tables.json and produces:
  1. NLDM table printout (for all 13 cells)
  2. Delay vs. Load plot for inverter family at tin=0.1225ns
  3. RC Model vs. SPICE comparison at midpoint (tin=0.1225ns, cload=0.0094pF)
  4. Saves plots as PDF/PNG
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# ─── Constants ──────────────────────────────────────────────
INPUT_TRANSITIONS_NS = [0.0100, 0.0231, 0.0531, 0.1225, 0.2823, 0.6507, 1.5000]
LOAD_CAPS_PF         = [0.0005, 0.0013, 0.0035, 0.0094, 0.0249, 0.0662, 0.1758]

MID_TIN_IDX   = 3   # 0.1225 ns
MID_CLOAD_IDX = 3   # 0.0094 pF
MID_TIN       = INPUT_TRANSITIONS_NS[MID_TIN_IDX]
MID_CLOAD     = LOAD_CAPS_PF[MID_CLOAD_IDX]

VDD  = 1.8
TEMP = 25

# ─── SKY130 TT parameters from Tutorial 2 ───────────────────
# nMOS W=0.42um, L=0.15um (smallest / x1 drive)
# Ref: sky130 PDK process specification
Kn       = 267e-6    # A/V² (process transconductance parameter for nMOS)
Vtn      = 0.48      # V    (nMOS threshold voltage)
Cox      = 8.46e-3   # F/m² (gate oxide capacitance density)
W_n      = 0.42e-6   # m
L_n      = 0.15e-6   # m
mu_n     = Kn / Cox  # carrier mobility (m²/Vs)
Ron_n    = 1.0 / (Kn * (W_n / L_n) * (VDD - Vtn))  # Ω

# Gate capacitance of smallest nMOS (Cg = Cox * W * L)
Cg_n     = Cox * W_n * L_n   # F

print(f"RC Model parameters:")
print(f"  Ron (nMOS x1) = {Ron_n:.1f} Ω")
print(f"  Cg  (nMOS x1) = {Cg_n*1e15:.3f} fF")

Ron_p_factor = 2.0   # PMOS is ~2× slower per unit width in SKY130
Ron_p        = Ron_n * Ron_p_factor


def rc_delay_fall(cload_pf):
    """RC model fall delay (ns) for invx1 at VDD=1.8V."""
    C = cload_pf * 1e-12
    return 0.69 * Ron_n * C * 1e9   # seconds -> ns


def rc_delay_rise(cload_pf):
    """RC model rise delay (ns) for invx1 at VDD=1.8V."""
    C = cload_pf * 1e-12
    return 0.69 * Ron_p * C * 1e9


def print_nldm_table(cell_name, cell_data, table_name):
    arr = np.array(cell_data["tables"][table_name])
    print(f"\n  {table_name}  [{cell_name}]")
    header = "  tin\\cload(pF) | " + " | ".join(f"{c:.4f}" for c in LOAD_CAPS_PF)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for i, tin in enumerate(INPUT_TRANSITIONS_NS):
        row_str = " | ".join(
            f"{arr[i,j]:.4f}" if not np.isnan(arr[i,j]) else "  NaN " 
            for j in range(len(LOAD_CAPS_PF))
        )
        print(f"  {tin:.4f} ns    | {row_str}")


def plot_inv_delay_vs_load(all_results, out_path):
    """Plot cell_rise and cell_fall delay vs load for inverter family at tin=0.1225ns."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Inverter Family: Delay vs. Load  (t_in = {MID_TIN} ns)", fontsize=13)

    inv_cells  = ["invx1", "invx2", "invx4", "invx8"]
    colors     = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    markers    = ["o", "s", "^", "D"]
    x          = LOAD_CAPS_PF

    for ax, table_name, title in zip(
        axes,
        ["cell_rise", "cell_fall"],
        ["Cell Rise Delay", "Cell Fall Delay"]
    ):
        for cell, color, marker in zip(inv_cells, colors, markers):
            if cell not in all_results:
                continue
            arr = np.array(all_results[cell]["tables"][table_name])
            y = arr[MID_TIN_IDX, :]
            ax.plot(x, y, marker=marker, color=color, label=cell, linewidth=1.8, markersize=6)

        ax.set_xscale("log")
        ax.set_xlabel("Load Capacitance (pF)", fontsize=11)
        ax.set_ylabel("Delay (ns)", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"\n✓ Plot saved: {out_path}")
    plt.close()


def rc_comparison(all_results, out_path):
    """
    Compare RC model vs SPICE for invx1 at all 7 load caps (tin fixed at midpoint).
    Also compute % error at the exact midpoint.
    """
    if "invx1" not in all_results:
        print("WARNING: invx1 data not found, skipping RC comparison.")
        return

    arr_rise = np.array(all_results["invx1"]["tables"]["cell_rise"])
    arr_fall = np.array(all_results["invx1"]["tables"]["cell_fall"])

    spice_rise = arr_rise[MID_TIN_IDX, :]
    spice_fall = arr_fall[MID_TIN_IDX, :]
    rc_rise    = np.array([rc_delay_rise(c) for c in LOAD_CAPS_PF])
    rc_fall    = np.array([rc_delay_fall(c) for c in LOAD_CAPS_PF])

    # ── Midpoint comparison ──────────────────────────────────
    mid_spice_rise = spice_rise[MID_CLOAD_IDX]
    mid_spice_fall = spice_fall[MID_CLOAD_IDX]
    mid_rc_rise    = rc_rise[MID_CLOAD_IDX]
    mid_rc_fall    = rc_fall[MID_CLOAD_IDX]

    print("\n" + "="*60)
    print("RC Model vs. SPICE Comparison")
    print(f"  Condition: tin = {MID_TIN} ns,  Cload = {MID_CLOAD} pF")
    print("="*60)
    print(f"  {'Metric':<22} {'RC Model (ns)':>14} {'SPICE (ns)':>12} {'Error %':>10}")
    print("  " + "-"*60)
    for label, rc_val, sp_val in [
        ("cell_rise (invx1)", mid_rc_rise, mid_spice_rise),
        ("cell_fall (invx1)", mid_rc_fall, mid_spice_fall),
    ]:
        if sp_val and not np.isnan(sp_val):
            err = (rc_val - sp_val) / sp_val * 100
            print(f"  {label:<22} {rc_val:>14.4f} {sp_val:>12.4f} {err:>+10.1f}%")
        else:
            print(f"  {label:<22} {rc_val:>14.4f}      (no SPICE data)")

    # ── Plot ────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"RC Model vs. SPICE — invx1  (t_in = {MID_TIN} ns)", fontsize=13)

    for ax, spice_y, rc_y, title in zip(
        axes,
        [spice_rise, spice_fall],
        [rc_rise,    rc_fall],
        ["Rise Delay", "Fall Delay"]
    ):
        ax.plot(LOAD_CAPS_PF, spice_y, "o-", color="#1f77b4", label="SPICE (NLDM)", linewidth=2)
        ax.plot(LOAD_CAPS_PF, rc_y,    "s--", color="#d62728", label="RC Model", linewidth=2)
        # Mark midpoint
        ax.axvline(MID_CLOAD, color="gray", linestyle=":", alpha=0.7, label=f"Cload={MID_CLOAD}pF")
        ax.set_xscale("log")
        ax.set_xlabel("Load Capacitance (pF)", fontsize=11)
        ax.set_ylabel("Delay (ns)", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.legend(fontsize=10)
        ax.grid(True, which="both", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"✓ RC comparison plot saved: {out_path}")
    plt.close()


def generate_markdown(all_results, out_path):
    """Write a markdown summary of all NLDM tables."""
    lines = ["# NLDM Characterization Report", "",
             "**Process:** SKY130  |  **Corner:** TT  |  **VDD:** 1.8V  |  **Temp:** 25°C", "",
             "## Index Vectors", "",
             "**Input Transition (ns):** " + str(INPUT_TRANSITIONS_NS), "",
             "**Load Capacitance (pF):** " + str(LOAD_CAPS_PF), ""]

    table_names = ["cell_rise", "cell_fall", "rise_transition", "fall_transition"]

    for cell_name in sorted(all_results.keys()):
        lines.append(f"## {cell_name}")
        for tname in table_names:
            arr = np.array(all_results[cell_name]["tables"][tname])
            lines.append(f"\n### {tname}")
            # Markdown table header
            header = "| t_in \\ C_load |" + "|".join(f" {c:.4f} " for c in LOAD_CAPS_PF) + "|"
            sep    = "|" + "|".join(["---"] * (len(LOAD_CAPS_PF) + 1)) + "|"
            lines.append(header)
            lines.append(sep)
            for i, tin in enumerate(INPUT_TRANSITIONS_NS):
                row = f"| {tin:.4f} ns |" + "|".join(
                    f" {arr[i,j]:.4f} " if not np.isnan(arr[i,j]) else " N/A "
                    for j in range(len(LOAD_CAPS_PF))
                ) + "|"
                lines.append(row)
        lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    print(f"✓ Markdown report saved: {out_path}")


def main():
    data_path = Path("results/nldm_tables.json")
    if not data_path.exists():
        print(f"ERROR: {data_path} not found. Run characterize.py first.")
        sys.exit(1)

    with open(data_path) as f:
        all_results = json.load(f)

    print(f"Loaded data for {len(all_results)} cells: {list(all_results.keys())}")
    out_dir = Path("results")

    # 1. Print all NLDM tables to console
    for cell_name, cell_data in all_results.items():
        print(f"\n{'='*70}")
        print(f"Cell: {cell_name}")
        for tname in ["cell_rise", "cell_fall", "rise_transition", "fall_transition"]:
            print_nldm_table(cell_name, cell_data, tname)

    # 2. Delay vs Load plot for inverter family
    plot_inv_delay_vs_load(all_results, out_dir / "inv_delay_vs_load.pdf")

    # 3. RC Model vs SPICE comparison
    rc_comparison(all_results, out_dir / "rc_vs_spice.pdf")

    # 4. Markdown report
    generate_markdown(all_results, out_dir / "nldm_report.md")


if __name__ == "__main__":
    main()
