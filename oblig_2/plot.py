import numpy as np
import matplotlib.pyplot as plt

# Load the data from the C++ simulation output file
sim_data = np.loadtxt("1d.txt")
r_sim = sim_data[:, 0].astype(int)
C_sim_complex = sim_data[:, 1] + 1j * sim_data[:, 2]
C_sim_real = np.real(C_sim_complex)

def C_analytic(T, r, L):
    '''Analytic expression of Correlation function'''
    J = -1
    beta = -1 / T
    exp_beta_J = np.exp(beta * J)
    term1_part1 = (1/3) * (exp_beta_J + 2)**r - (1/3) * (exp_beta_J - 1)**r
    term1_part2 = (1/3) * (exp_beta_J + 2)**(L - r) - (1/3) * (exp_beta_J - 1)**(L - r)
    term1 = 3 * np.exp((2j * np.pi) / 3) * term1_part1 * term1_part2

    term2_part1 = (1/3) * (exp_beta_J + 2)**r - (1/3) * (exp_beta_J - 1)**r
    term2_part2 = (1/3) * (exp_beta_J + 2)**(L - r) - (1/3) * (exp_beta_J - 1)**(L - r)
    term2 = 3 * np.exp((-2j * np.pi) / 3) * term2_part1 * term2_part2

    term3_part1 = (2/3) * (exp_beta_J - 1)**r + (1/3) * (exp_beta_J + 2)**r
    term3_part2 = (2/3) * (exp_beta_J - 1)**(L - r) + (1/3) * (exp_beta_J + 2)**(L - r)
    term3 = 3 * term3_part1 * term3_part2

    denominator = (2 * (exp_beta_J - 1)**L + (exp_beta_J + 2)**L)
    return (term1 + term2 + term3) / denominator

L = 16
T = 0.25 # Has to be manually set for each T
r_analytic = np.arange(0, L) 
C_exact = C_analytic(T, r_analytic, L)

# Plotting
plt.plot(r_sim, C_sim_real, "o", label="Simulation", color='b')
plt.plot(r_analytic, np.real(C_exact), label="Exact", color='r', linestyle='-')
plt.xlabel('r')
plt.ylabel('Real part of C(r)')
plt.legend()
plt.title(f'Real part of Correlation C(r) at T={T}, L={L}')
plt.grid(True)
plt.savefig(f"wolff_1D_T_{T:.2f}_real.pdf")
plt.show()