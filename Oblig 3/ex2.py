import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

sigma_T = 6.65e-25
omega_lambda = 0.692
omega_m = 0.308
omega_r = 0

c = 2.99e10

Mpc = 1e3*3.086e18

H_0 = 72e3 /Mpc

n_e  = lambda z :  1.9e-7*(1+z)**3
H = lambda z : H_0*np.sqrt(omega_lambda + omega_m*(1+z)**3)


integrand = lambda z : c*n_e(z)*sigma_T/(H(z)*(1+z))


N = 1000
dz = 10/N
z = np.linspace(0,10,N)
tau = np.zeros(N)


for i in range(1,N):
    tau[i] = tau[i-1] + integrand(z[i-1])*dz

plt.plot(z, tau)
plt.title(r"Optical Depth $\tau_e$")
plt.xlabel(r"Redshift $z$")
plt.ylabel(r"$\tau_e$")
plt.show()
