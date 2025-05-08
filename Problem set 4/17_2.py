import matplotlib.pyplot as plt
import numpy as np

V = np.linspace(1,10,100)

a = 1
b = 2
N = 1
k_b = 1

P_c = a/(27*b**2)
T_c = 8*a/(27*k_b*b)
print(T_c)
V_c = 3*N*b

T = T_c

P_vdW =  N*k_b*T/(V-b*N) - a*N**2/V**2
P_ideal = N*k_b*T/V

plt.figure(figsize=(7,5))
plt.plot(V/V_c,P_vdW/P_c, label='Van der Waal gas')
plt.plot(V/V_c,P_ideal/P_c, label = 'Ideal gas')
plt.xlabel("V")
plt.ylabel("P")
plt.title(f'Pressure over V at: T={T}, N={N}, a={a}, b={b}.')
plt.grid()
plt.legend()
plt.show()

# T = np.linspace(0,10,100)
# Ps = np.array([P_c-1,P_c,P_c+1])
# plt.figure(figsize=(7,5))
# for P in Ps:
#     T = P*(V-b*N)/(N*k_b) + a*N*N*(V-b*N)/(V*V*N*k_b)
#     plt.plot(V,T,label=f'P:{P}')

# plt.title(f'T over V at: N={N}, a={a}, b={b}.')
# plt.xlabel("V")
# plt.ylabel("T")
# plt.grid()
# plt.legend()
# plt.show()



# V = np.linspace(0,10,100)
# P_vdW =  N*k_b*T/(V-b*N) - a*N**2/V**2
# mu = -k_b*T*np.log(V-b*N) + N*k_b*T*b/(V-b*N) + k_b*T*np.log(N) + k_b*T - 2*a*N/V

# plt.figure(figsize=(7,5))
# plt.plot(P_vdW,mu)
# plt.title(f'mu over P.')
# plt.xlabel("P")
# plt.ylabel("mu")
# plt.grid()
# plt.show()
