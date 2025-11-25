import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
order = "o4"
g="g1"
kind = "SE"       
tot = 21
all_data = []
beta =10
U = 1
# Load all data into a list
for i in [0,1,2]:
# for i in [1,13,14]:
    # fname = f"data/Summed_data/{i}i_shot_{kind}.txt"
    fname = f"det/data/Individual_data/SE_{order}_{g}_n{i}_i1.txt"
    # fname = "data/Summed_data/Summed_sigma_ac.txt"
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
qx_val = 0
qy_val = np.pi

# Use np.isclose to filter
# df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) &  (df["wn"] >= 0) ]
df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) ]
# df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) ]
# ==== plot ====


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)


for i in sorted(df_filtered["iter"].unique()):
    sub = df_filtered[df_filtered["iter"] == i]
    ax1.plot(sub["wn"], sub["Re"],".-", label=f"SE_{order}_{g}_n{i}")
    ax2.plot(sub["wn"], sub["Im"], ".-", label=f"SE_{order}_{g}_n{i}")

        


ax1.set_ylabel(" Re $\Sigma$ ",fontsize=11)
ax2.set_ylabel("Im $\Sigma$ ",fontsize=11)
ax2.set_xlabel(" DLR Matsubara Frequency $\\omega_n$",fontsize=11)


ax1.legend(ncol=3)
ax2.legend(ncol=3)
ax1.grid(True)
ax2.grid(True)
plt.suptitle(r'GF2, $q=(\pi,\pi), N_c=4,\beta=4, U=1$',fontsize=13)
plt.tight_layout()
plt.show()


last_iter =1

# Select only that iteration
sub_last = df_filtered[df_filtered["iter"] == last_iter]

# Stack wn and Im into two columns
out = np.column_stack((sub_last["wn"], sub_last["Im"]))

# # Save to disk
# np.savetxt(f"data/o3_DCA_N0_NC16.txt",
#            out,
#            fmt="%.8e")


