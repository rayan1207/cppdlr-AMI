#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-figure, two-panel scaling plots for cppdlr+AMI
Left: Order-4 FK & CK strong-scaling vs CPUs with ideal 1/p
Right: Order-3 CK vs L at 64 CPUs
"""

import numpy as np
import matplotlib.pyplot as plt

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

# ---------------- Data (as provided) ----------------
# ---- Order 4 (L = 4) strong-scaling vs CPUs ----
cpu_o4_L4     = np.array([16, 32, 64, 96, 128, 192])
Fk_fixed_L4o4 = np.array([5123165, 2395940, 1183716, 939290, 694375, 444012])  # ms
CT_fixed_L4o4 = np.array([  83344,   39985,    19938,   14721,   11664,   7559])  # ms

# ---- Order 3 (CPU = 64) CT vs L ----
L_o3               = np.array([2, 4, 6, 8, 10, 12])
CT_fixed_cpu64_o3  = np.array([31, 4666, 76408, 645652, 3666784, 13733880])  # ms

# ---------------- Helpers ----------------
def ideal_time_series(cpus, t0, p0):
    """Ideal strong scaling: T(p) = T(p0) * (p0 / p)"""
    cpus = cpus.astype(float)
    return t0 * (p0 / cpus)

# Reference CPU count (use the smallest available as baseline)
p0 = cpu_o4_L4[0]

# Ideal lines (match each series' baseline)
fk_ideal = ideal_time_series(cpu_o4_L4, Fk_fixed_L4o4[0], p0)
ck_ideal = ideal_time_series(cpu_o4_L4, CT_fixed_L4o4[0], p0)

# ---------------- Plot ----------------
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

# --- Left panel: Order 4 strong scaling ---
ax = axes[0]
ax.set_title(r'Order 4 (L = 2): $CK$ Wall-time scaling ')
# Measured
ax.plot(cpu_o4_L4, Fk_fixed_L4o4, marker="o", linewidth=1.8, markersize=5,
        label=r'$K$ (measured)', color="tab:blue")
ax.plot(cpu_o4_L4, CT_fixed_L4o4, marker="s", linewidth=1.8, markersize=5,
is         label=r'$C$ (measured)', color="tab:orange")
# Ideals
ax.plot(cpu_o4_L4, fk_ideal, linestyle="--", linewidth=1.2, label=r' $K$ ideal (1/p)', color="tab:blue", alpha=0.6)
ax.plot(cpu_o4_L4, ck_ideal, linestyle="--", linewidth=1.2, label=r' $C$ ideal (1/p)', color="tab:orange", alpha=0.6)
''
ax.set_xscale("log", base=2)
ax.set_yscale("log")
ax.set_xlabel("CPUs")
ax.set_ylabel("Time (ms)")
# ax.grid(True, which="both")
ax.legend(loc="best", frameon=False)

# --- Right panel: Order 3 CT vs L @ 64 CPUs ---
ax2 = axes[1]
ax2.set_title("Order 3: C vs L @ 64 CPUs")
ax2.plot(L_o3, CT_fixed_cpu64_o3, marker="^", linewidth=1.8, markersize=5,
         label="C (measured)", color="tab:orange")
ax2.set_yscale("log")
ax2.set_xlabel("L")
ax2.set_ylabel("Time (ms)")
# ax2.grid(True, which="both")
ax2.legend(loc="best", frameon=False)

# Overall figure title (optional)
# fig.suptitle("cppdlr+AMI Scaling", fontsize=13, fontweight="bold")

# Save outputs
# plt.savefig("cppdlr_AMI_scaling_panels.png")
plt.savefig("cppdlr_AMI_scaling_panels.pdf")
plt.show()

print("✅ Saved: cppdlr_AMI_scaling_panels.png and cppdlr_AMI_scaling_panels.pdf")
