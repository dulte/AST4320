import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
try:
    import seaborn as sns
    sns.set_style("darkgrid")
except:
    print("Everything Looks Nicer With Seaborn!")


def T_photon(a):
    return 4*1./a
def T_gas(a):
    T = np.where(a<1./(1100+1),T_photon(a),4e-3*1./a**(2))
    return T


def jeans_length_before(z):
    c = constants.c
    G = constants.G
    H0 = 2e-18 #per sec
    rho_0 = 3*H0**2/(8*np.pi*G)

    return c/(np.sqrt(3))*np.sqrt(np.pi/(G*rho_0))*(1+z)**(-3./2)


def jeans_length_after(z):
    G = constants.G
    H0 = 2e-18 #per sec
    rho_0 = 3*H0**2/(8*np.pi*G)
    mp = constants.m_p
    k = constants.k
    mu = 0.62

    return 4*10**(-3/2)*(1+z)**(-1./2)*np.sqrt(k/(mu*mp))*np.sqrt(np.pi/(G*rho_0))

def jeans_mass_before(z):
    c = constants.c
    G = constants.G
    H0 = 2e-18 #per sec
    rho_0 = 3*H0**2/(8*np.pi*G)

    return (np.pi**(5./2))/(6*G**(3/2.)*rho_0**(1./2))*(c/np.sqrt(3))**3*(1+z)**(-3./2)

def jeans_mass_after(z):
    G = constants.G
    H0 = 2e-18 #per sec
    rho_0 = 3*H0**2/(8*np.pi*G)

    return (np.pi**(5./2))/(6*G**(3/2.)*rho_0**(1./2))*4**3*10**(0.5)*(1+z)**(3./2)

def jeans_length(z):
    return np.where(z>1100,jeans_length_before(z),jeans_length_after(z))


def jeans_mass(z):
    return np.where(z>1100,jeans_mass_before(z),jeans_mass_after(z))



if __name__ == '__main__':
    #a)

    a = np.linspace(1e-4,1,1000)
    plt.loglog(a,T_photon(a),label="Photons")
    plt.loglog(a,T_gas(a),label="Gas")
    plt.title("Temperature of gas and photons",fontsize=20)
    plt.xlabel("a",fontsize=15)
    plt.ylabel("T [K]",fontsize=15)
    plt.legend(loc="best",fontsize=15)
    plt.show()

    #b)
    z = np.linspace(2000,0,10000)
    jeans_length_gas = jeans_length(z)
    jeans_mass = jeans_mass(z)
    plt.loglog(z,jeans_length_gas)
    plt.xlim(z[0],z[-1])
    plt.title("Jeans Length Before and After Decoupling",fontsize=20)
    plt.xlabel("z",fontsize=15)
    plt.ylabel(r"$\lambda_J$",fontsize=15)
    plt.show()


    plt.loglog(z,jeans_mass)
    plt.title("Jeans Mass Before and After Decoupling",fontsize=20)
    plt.xlabel("z",fontsize=15)
    plt.ylabel(r"$M_J$",fontsize=15)
    plt.xlim(z[0],z[-1])
    plt.show()
