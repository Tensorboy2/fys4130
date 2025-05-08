import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

L_values = [8, 16, 32]
colors = ['r', 'g', 'b']
markers = ['o', 's', '^']

fig, ax = plt.subplots(figsize=(8, 6))

for i, L in enumerate(L_values):
    filename = f"binder_gamma_L{L}_q3.dat"
    try:
        data = np.loadtxt(filename)
        T_over_J = data[:, 0]
        Gamma = data[:, 1]
        ax.plot(T_over_J, Gamma, marker='.', linestyle='-', color=colors[i], label=f"L={L}")
    except FileNotFoundError:
        print(f"File not found: {filename}")

# Vertical line at Tc
Tc = 0.995
ax.axvline(x=Tc, ymin=0, ymax=1, color="black", linestyle='--', label="$T_c=0.9950$")

# Main plot labels
ax.set_xlabel("T/J")
ax.set_ylabel(r"$\Gamma$")
ax.set_title(r"$\Gamma$ vs. Temperature for different $L$")
ax.legend()
ax.grid(True)

# Add zoomed-in inset
axins = inset_axes(ax, width="40%", height="40%", loc="upper right")
zoom_xlim = (0.93, 1.02)
zoom_ylim = (1.0, 1.4)

for i, L in enumerate(L_values):
    filename = f"binder_gamma_L{L}_q3.dat"
    try:
        data = np.loadtxt(filename)
        T_over_J = data[:, 0]
        Gamma = data[:, 1]
        axins.plot(T_over_J, Gamma, marker='.', linestyle='-', color=colors[i])
    except FileNotFoundError:
        continue

axins.set_xlim(zoom_xlim)
axins.set_ylim(zoom_ylim)
axins.axvline(x=Tc, color="black", linestyle='--')
axins.set_xticks(np.round(np.linspace(*zoom_xlim, 3), 3))
axins.set_yticks(np.round(np.linspace(*zoom_ylim, 3), 2))
axins.grid(True)

# Save and show
plt.savefig("binder_cumulant_comparison_zoomed.pdf")
plt.show()
