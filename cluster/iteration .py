#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kind = "SE"       
tot = 21
all_data = []
beta =10
U = 2
ilist = [0,1,2,3,19]
# Load all data into a list
for i in ilist:
# for i in [1,13,14]:
    fname = f"GF2/U{U}/data/Summed_data/{i}i_shot_{kind}.txt"
    # fname = "DCA/data/Summed_data/single_shot_AC_ph.txt"
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
df_filtered = df[np.isclose(df["qx"], 0) & np.isclose(df["qy"], np.pi)&  (df["wn"] >= 0)  & (df["wn"] < 20)  ]
df_filtered1 = df[np.isclose(df["qx"], np.pi) & np.isclose(df["qy"], np.pi)&  (df["wn"] >= 0) & (df["wn"] < 20) ]
df_filtered2 = df[np.isclose(df["qx"], 0) & np.isclose(df["qy"], 0)&  (df["wn"] >= 0) & (df["wn"] < 20) ]
# df_filtered = df[np.isclose(df["qx"], qx_val) & np.isclose(df["qy"], qy_val) ]
# ==== plot ====

plt.figure(figsize=[15,7])
plt.suptitle(f'DCA-IPT4 solver, U={U}',fontsize=18)
for i in sorted(df_filtered["iter"].unique()):
     sub = df_filtered[df_filtered["iter"] == i]
     sub1 = df_filtered1[df_filtered1["iter"] == i]
     sub2 = df_filtered2[df_filtered2["iter"] == i]
     plt.subplot(1,3,1)
     plt.title(r'$k=[0,\pi]$',fontsize=15)
     plt.plot(sub["wn"], sub["Im"], "^-", label=f"i={i}")
     plt.ylabel(r'Im $\Sigma$',fontsize=12)
     plt.xlabel(r'$i\omega$',fontsize=12)
     plt.legend(ncol=3)
     plt.grid(True)
     # plt.ylim(-8,-0.4)
     plt.subplot(1,3,2)
     plt.title(r'$k=[\pi,\pi]$',fontsize=15)
     plt.plot(sub1["wn"], sub1["Im"], "^-", label=f"i={i}")
     plt.grid(True)
     plt.ylabel(r'$i\omega$',fontsize=12)
     plt.subplot(1,3,3)
     plt.title(r'$k=[0,0]$',fontsize=15)
     plt.plot(sub2["wn"], sub2["Im"], "^-", label=f"i={i}")
     plt.grid(True)
     plt.ylabel(r'$i\omega$',fontsize=12)
plt.tight_layout()
# plt.savefig(f'IPT4_U{U}_Im.pdf',dpi=400)
plt.show()


# plt.suptitle('DCA-IPT4 solver, U=2',fontsize=18)
plt.figure(figsize=[15,7])
for i in sorted(df_filtered["iter"].unique()):
     sub = df_filtered[df_filtered["iter"] == i]
     sub1 = df_filtered1[df_filtered1["iter"] == i]
     sub2 = df_filtered2[df_filtered2["iter"] == i]
     plt.subplot(1,3,1)
     plt.title(r'$k=[0,\pi]$',fontsize=15)
     plt.plot(sub["wn"], sub["Re"], "^-", label=f"i={i}")
     plt.ylabel(r'RE $\Sigma$',fontsize=12)
     plt.xlabel(r'$i\omega$',fontsize=12)
     plt.legend(ncol=3)
     plt.grid(True)
     # plt.ylim(-1e-04,1e-04)
     plt.subplot(1,3,2)
     plt.title(r'$k=[\pi,\pi]$',fontsize=15)
     plt.plot(sub1["wn"], sub1["Re"], "^-", label=f"i={i}")
     plt.grid(True)
     plt.ylabel(r'$i\omega$',fontsize=12)
     plt.subplot(1,3,3)
     plt.title(r'$k=[0,0]$',fontsize=15)
     plt.plot(sub1["wn"], sub2["Re"], "^-", label=f"i={i}")
     plt.grid(True)
     plt.ylabel(r'$i\omega$',fontsize=12)
plt.tight_layout()
# plt.savefig(f'IPT4_U{U}_RE.pdf',dpi=400)
plt.show()
# plt.close()
     
     
     
#    
    



# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)


# for i in sorted(df_filtered["iter"].unique()):
#     sub = df_filtered[df_filtered["iter"] == i]
#     ax1.plot(sub["wn"], sub["Re"],".-", label=f"i={i}")
#     ax2.plot(sub["wn"], sub["Im"]/(U**2), ".-", label=f"i={i}")
#     # ax2.set_ylim(-0.010,0)

        


# ax1.set_ylabel(" Re $\Sigma$ ",fontsize=11)
# ax2.set_ylabel("Im $\Sigma$ ",fontsize=11)
# ax2.set_xlabel(" DLR Matsubara Frequency $\\omega_n$",fontsize=11)


# ax1.legend(ncol=3)
# ax2.legend(ncol=3)
# ax1.grid(True)
# ax2.grid(True)
# plt.suptitle(r'GF2, $q=(\pi,0), N_c=64,\beta=33, U=13$',fontsize=13)
# plt.tight_layout()
# plt.show()


last_iter =ilist[-1]

# Select only that iteration
sub_last = df_filtered[df_filtered["iter"] == last_iter]

# Stack wn and Im into two columns
out = np.column_stack((sub_last["wn"], sub_last["Im"]))

print(last_iter)

# plt.plot(sub_last["wn"], sub_last["Im"],'^-')
# plt.show()
# Save to disk
# np.savetxt(f"data/o4_DCA_Nc_n1_beta5_U{U}.txt",
#             out,
#             fmt="%.8e")


