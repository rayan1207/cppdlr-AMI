import numpy as np
import matplotlib.pyplot as plt

# =========================
# Load data
# =========================
# 4th order
o4_u2   = np.loadtxt('o4_DCA_Nc_n1_beta5_U2.txt')
o4_u3   = np.loadtxt('o4_DCA_Nc_n1_beta5_U3.txt')
o4_u4   = np.loadtxt('o4_DCA_Nc_n1_beta5_U4.txt')
o4_u45  = np.loadtxt('o4_DCA_Nc_n1_beta5_U4.5.txt')
o4_u5   = np.loadtxt('o4_DCA_Nc_n1_beta5_U5.txt')
o4_u525 = np.loadtxt('o4_DCA_Nc_n1_beta5_U5.25.txt')

# 2nd order
# NOTE: change the filenames here if your actual 2nd-order files are named differently
o2_u2   = np.loadtxt('o2_DCA_Nc_n1_beta5_U2.txt')
o2_u3   = np.loadtxt('o2_DCA_Nc_n1_beta5_U3.txt')
o2_u4   = np.loadtxt('o2_DCA_Nc_n1_beta5_U4.txt')
o2_u45   = np.loadtxt('o2_DCA_Nc_n1_beta5_U4.5.txt')
o2_u5   = np.loadtxt('o2_DCA_Nc_n1_beta5_U5.txt')
o2_u6   = np.loadtxt('o2_DCA_Nc_n1_beta5_U6.txt')
o2_u8   = np.loadtxt('o2_DCA_Nc_n1_beta5_U8.txt') 
o2_u10   = np.loadtxt('o2_DCA_Nc_n1_beta5_U10.txt')   # use U=8 or 8.25 as you like
# o2_u10  = np.loadtxt('o2_DCA_Nc_n1_beta5_U10.txt')   # loaded if you need later

# =========================
# Extract columns
# =========================
# 4th order
wn_o4_u2,   im_o4_u2   = o4_u2[:, 0], o4_u2[:, 1]
wn_o4_u3,   im_o4_u3   = o4_u3[:, 0], o4_u3[:, 1]
wn_o4_u4,   im_o4_u4   = o4_u4[:, 0], o4_u4[:, 1]
wn_o4_u45,  im_o4_u45  = o4_u45[:, 0], o4_u45[:, 1]
wn_o4_u5,   im_o4_u5   = o4_u5[:, 0], o4_u5[:, 1]
wn_o4_u525, im_o4_u525 = o4_u525[:, 0], o4_u525[:, 1]

# 2nd order
wn_o2_u2,  im_o2_u2  = o2_u2[:, 0], o2_u2[:, 1]
wn_o2_u3,  im_o2_u3  = o2_u3[:, 0], o2_u3[:, 1]
wn_o2_u4,  im_o2_u4  = o2_u4[:, 0], o2_u4[:, 1]
wn_o2_u45,  im_o2_u45  = o2_u45[:, 0], o2_u45[:, 1]
wn_o2_u5,  im_o2_u5  = o2_u5[:, 0], o2_u5[:, 1]
wn_o2_u6,  im_o2_u6  = o2_u6[:, 0], o2_u6[:, 1]
wn_o2_u8,  im_o2_u8  = o2_u8[:, 0], o2_u8[:, 1]
wn_o2_u10, im_o2_u10 = o2_u10[:, 0], o2_u10[:, 1]

# =========================
# Matplotlib look
# =========================
plt.rcParams['figure.dpi'] = 120
plt.rcParams['font.size'] = 12

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),)

# =========================
# First subplot:
# o4: U = 2, 3, 4, 4.5 (solid markers)
# o2: U = 2, 3, 4, 5   (open markers)
# =========================

# 4th order (solid markers)
ax1.plot(wn_o4_u2,  im_o4_u2,  '-o',  label=r'$U=2$ (4th)',  markersize=4)
ax1.plot(wn_o4_u3,  im_o4_u3,  '-s',  label=r'$U=3$ (4th)',  markersize=4)
ax1.plot(wn_o4_u4,  im_o4_u4,  '-^',  label=r'$U=4$ (4th)',  markersize=4)
ax1.plot(wn_o4_u45, im_o4_u45, '-D',  label=r'$U=4.5$ (4th)', markersize=4)

# 2nd order (open markers: markerfacecolor='none')
ax1.plot(wn_o2_u2,  im_o2_u2,  '--o', label=r'$U=2$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax1.plot(wn_o2_u3,  im_o2_u3,  '--s', label=r'$U=3$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax1.plot(wn_o2_u4,  im_o2_u4,  '--^', label=r'$U=4$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax1.plot(wn_o2_u4,  im_o2_u45,  '--^', label=r'$U=4.5$ (2nd)',  markersize=4,
         markerfacecolor='none')


ax1.set_xlabel(r'$\omega_n$')
ax1.set_ylabel(r'$\mathrm{Im}\,\Sigma(i\omega_n)$')
ax1.set_title('Moderate $U$ comparison')
ax1.grid(True, alpha=0.3)
ax1.legend(frameon=False)

# =========================
# Second subplot:
# Metal/insulator:
#   o4: U = 4.5, 5, 5.25
#   o2: U = 6, 8 (here 8.25 in file)
# =========================

# 4th order (solid markers)
ax2.plot(wn_o4_u45,  im_o4_u45,  '-o', label=r'$U=4.5$ (4th)', markersize=4)
ax2.plot(wn_o4_u5,   im_o4_u5,   '-s', label=r'$U=5$ (4th)',   markersize=4)
ax2.plot(wn_o4_u525, im_o4_u525, '-^', label=r'$U=5.25$ (4th)', markersize=4)

# 2nd order (open markers)
ax2.plot(wn_o2_u5,  im_o2_u5,  '--D', label=r'$U=5$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax2.plot(wn_o2_u6,  im_o2_u6,  '--D', label=r'$U=6$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax2.plot(wn_o2_u8,  im_o2_u8,  '--v', label=r'$U=8$ (2nd)',  markersize=4,
         markerfacecolor='none')
ax2.plot(wn_o2_u10,  im_o2_u10,  '--v', label=r'$U=10$ (2nd)',  markersize=4,
         markerfacecolor='none')


ax2.set_xlabel(r'$\omega_n$')
ax2.set_title('Metal–insulator behavior')
ax2.grid(True, alpha=0.3)
ax2.legend(frameon=False)

plt.tight_layout()
plt.savefig('IPT4_IPT2.pdf',dpi=800)
plt.show()
