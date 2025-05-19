import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kind = "SE"       
tot = 15        
all_data = []
beta =10
U = 4
# Load all data into a list

for i in [1,]:
    fname = f"{i}i_shot_{kind}.txt"
    try:
        data = np.loadtxt(fname, comments="#")
        for row in data:
            wn, qx, qy, re, im = row
            all_data.append([i, wn, qx, qy, re, im])
    except FileNotFoundError:
        print(f"Missing: {fname}")

# Convert to DataFrame
df = pd.DataFrame(all_data, columns=["iter", "wn", "qx", "qy", "Re", "Im"])

# ==== user input ====
qx_val = np.pi/2
qy_val = np.pi/2

# Use np.isclose to filter
df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val)]

# ==== plot ====
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

for i in sorted(df_filtered["iter"].unique()):
    sub = df_filtered[df_filtered["iter"] == i]
    ax1.plot(sub["wn"], sub["Re"],".-", label=f"i={i}")
    ax2.plot(sub["wn"], sub["Im"], ".-", label=f"i={i}")

ax1.set_ylabel("Re")
ax2.set_ylabel("Im")
ax2.set_xlabel(" DLR Matsubara Frequency $\\omega_n$")

ax1.set_title(f"{kind} at q ≈ ({qx_val:.3f}, {qy_val:.3f}), L=15x15,beta={beta},U={U}")
ax1.legend(ncol=3)
ax2.legend(ncol=3)
ax1.grid(True)
ax2.grid(True)
plt.tight_layout()
# plt.savefig('qpp_beta10_U4.pdf',dpi=500)
plt.show()
