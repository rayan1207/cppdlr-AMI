import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
order = "o4"
g="g0"
kind = "SE"       
tot = 9
all_data = []
beta =10
U = 1
trial=2
# Load all data into a list
for i in [0,1,2]:
# for i in [1,13,14]:
    # fname = f"data/Summed_data/{i}i_shot_{kind}.txt"
    fname = f"DCA/U1/data/Individual_data/SE_{order}_{g}_n{i}_i0.txt"
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
qx_val = np.pi
qy_val = np.pi

# Use np.isclose to filter
# df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) &  (df["wn"] >= 0) ]
df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) &    (df["wn"] <20 ) & (df["wn"] >-20)  ]
# df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) ]
# ==== plot ====


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
label_o4 =["SE_o4_g0_n2","SE_o4_g1_n2","SE_o4_g2_n0"]

indx=0
for i in sorted(df_filtered["iter"].unique()):
    sub = df_filtered[df_filtered["iter"] == i]
    ax1.plot(sub["wn"], sub["Re"],".-", label=label_o4[indx])
    ax2.plot(sub["wn"], sub["Im"], ".-", label=label_o4[indx])
    indx =indx+1

        


ax1.set_ylabel(" Re $\Sigma$ ",fontsize=11)
ax2.set_ylabel("Im $\Sigma$ ",fontsize=11)
ax2.set_xlabel(" DLR Matsubara Frequency $\\omega_n$",fontsize=11)


ax1.legend(ncol=2)
ax2.legend(ncol=2)
ax1.grid(True)
ax2.grid(True)
# plt.suptitle(f"single-shot {order}, trial-{trial} " + r'$\beta=5$, U=1, $k=[\pi,\pi]$ ',fontsize=13)
plt.suptitle(r'Three fourth order diagram considered $k=[\pi,\pi]$')
plt.tight_layout()
plt.savefig("Eff_3o4_SE_pp.pdf",dpi=300)
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


