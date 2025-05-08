import matplotlib.pyplot as plt
import numpy as np

# Plot for <m>
data_m = np.loadtxt("magnetization_L16_q3_2.dat")
plt.plot(data_m[:, 0], np.real(data_m[:, 1]), marker='o', linestyle='-')
plt.xlabel("T/J")
plt.ylabel("<m>")
plt.title("Average Magnetization per Site (L=16)")
plt.grid(True)
plt.savefig("wolff_2D_avg_m_over_T.pdf")
plt.show()

# Plot for <|m|^2>
data_m2 = np.loadtxt("magnetization_squared_L16_q3_2.dat")
plt.plot(data_m2[:, 0], np.real(data_m2[:, 1]), marker='o', linestyle='-')
plt.xlabel("T/J")
plt.ylabel("<|m|^2>")
plt.title("Average Magnetization Squared per Site (L=16)")
plt.grid(True)
plt.savefig("wolff_2D_avg_m_squared_over_T.pdf")
plt.show()