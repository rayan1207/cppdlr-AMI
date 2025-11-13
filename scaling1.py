#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Proposal-quality scaling plots for cppdlr+AMI
One figure per order (2, 3, 4), with FK and CK subplots
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext
# ---------------- Appearance ----------------
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.alpha": 0.35,
})

# ---------------- Data ----------------
# ---- Order 2 (L=4) ----
o2_Elist = 1331
cpu_o2_L4 = np.array([2, 8, 16, 24, 48, 64])
Fk_fixed_L4o2 = np.array([95, 25, 13.3, 10, 5.5, 4])               # ms
CT_fixed_L4o2 = np.array([41, 11, 5.3, 3, 1.7, 1.3])                  # ms (shorter)

# ---- Order 3 (L=4) ----
o3_Elist = 161051
cpu_o3_L4 = np.array([2, 8, 16, 24, 36, 48, 64])
Fk_fixed_L4o3 = np.array([51582, 12773, 6378, 4203, 3005, 2169, 1628])  # ms
CT_fixed_L4o3 = np.array([129272, 32200, 15971, 10664, 7136, 5342, 3366]) 

# ---- Order 4 (L=4) ----
o3_Elist = 17569200
cpu_o4_L4 = np.array([16, 32, 64, 96, 128, 192])
Fk_fixed_L4o4 = np.array([5123165, 2395940, 1183716, 939290, 694375, 444012]) # ms
CT_fixed_L4o4 = np.array([83344, 39985, 19938, 14721, 11664, 7559])           # ms


L_o3 = np.array([2, 4, 6, 8, 10, 12])
CT_fixed_cpu64_o3 = np.array([31, 4666, 76408, 645652, 3666784, 13733880])   # ms
L_o2 = np.array([4, 8, 12, 16, 24, 32])
CT_fixed_cpu64_o2 = np.array([1.3, 15, 520, 2635, 27571, 151474])           # ms

# ---------------- Helpers ----------------
def ideal_time_series(cpus, t0, p0):
    """Ideal strong scaling: T(p) = T(p0)*(p0/p)"""
    return t0 * (p0 / cpus.astype(float))

def plot_order(cpu, fk, ck, order):
    """Make one figure for each order with M(& CK subplots"""
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=False)
    fig.suptitle(f"Order {order} Scaling (L = 2)", fontsize=13, fontweight="bold")

    # --- FK subplot ---
    ax_fk = axes[0]
    ax_fk.plot(cpu, np.log10(fk), marker="o", linewidth=1.8, markersize=5, label="Measured FK", color="tab:blue")
    ax_fk.plot(cpu, np.log10(ideal_time_series(cpu, fk[0], cpu[0])),
               linestyle="--", linewidth=1.2, color="tab:gray", label="Ideal 1/p")
    ax_fk.set_xscale("log", base=2)
    # ax_fk.set_yscale("log")
    ax_fk.set_xlabel("CPUs")
    ax_fk.set_ylabel(" Log Time (ms)")
    ax_fk.set_title("Frequency Kernel (CK)")
    
    # # Simple log ticks - only show 10^n labels
    # ax_fk.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    # ax_fk.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    ax_fk.grid(True, which="both", alpha=0.4)
    ax_fk.legend()

    # --- CK subplot ---
    ax_ck = axes[1]
    n = min(len(cpu), len(ck))
    ax_ck.plot(cpu, np.log10(ck), marker="s", linewidth=1.8, markersize=5, label="Measured CK", color="tab:orange")
    ax_ck.plot(cpu, np.log10(ideal_time_series(cpu, ck[0], cpu[0])),
               linestyle="--", linewidth=1.2, color="tab:gray", label="Ideal 1/p")
    ax_ck.set_xscale("log", base=2)
    ax_ck.set_xlabel("CPUs")
    ax_ck.set_title("Coefficient Kernel (CT)")
    ax_ck.set_ylabel(" Log Time (ms)")
    ax_ck.grid(True, which="both", alpha=0.4)
    ax_ck.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"scaling_order_{order}_FK_CK.png")
    plt.show()

# ---------------- Plot for each order ----------------
# plot_order(cpu_o2_L4, Fk_fixed_L4o2, CT_fixed_L4o2, order=2)
# plot_order(cpu_o3_L4, Fk_fixed_L4o3, CT_fixed_L4o3, order=3)
plot_order(cpu_o4_L4, Fk_fixed_L4o4, CT_fixed_L4o4, order=4)

print("✅ Saved: scaling_order_2_FK_CK.png, scaling_order_3_FK_CK.png, scaling_order_4_FK_CK.png")

plt.figure(figsize=(9, 4.2))
# ---- Order 2 ----
plt.subplot(1, 2, 1)
plt.plot(L_o2, CT_fixed_cpu64_o2/1000, marker="o", linewidth=2.4, markersize=5,color='tab::blue')
plt.xlabel("L")
plt.ylabel("CK time (s)")
plt.title("Order 2: CT time vs L @ 64 CPUs")
plt.grid(True, which="both")

# ---- Order 3 ----
plt.subplot(1, 2, 2)
plt.plot(L_o3, CT_fixed_cpu64_o3/1000, marker="^", linewidth=2.4, markersize=5,color='tab::orange')
plt.xlabel("L")
plt.title("Order 3: CT time vs L @ 64 CPUs")
plt.grid(True, which="both")

plt.tight_layout()
plt.savefig("cpu64_CK_time_vs_L_combined.png", dpi=300)
plt.show()