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
    return 4e-3*1./a**(2)

def sound_speed(z,type_):
    c = constants.c
    k = constants.k
    mp = constants.m_p
    mu = 0.63
    if type_=="photons":
        return c/np.sqrt(3)
    else:
        return np.sqrt(T_gas(1./(1+z))*k/(mu*mp))

def density(z):
    G = constants.G
    a = 1./(z+1)
    H = 2e-18*np.sqrt(0.3*a**(-3)+0.7)
    rho = 3*H**2/(8*np.pi*G)
    return rho



def jeans_length(z,type_):
    G = constants.G
    c = np.zeros_like(z)
    c = np.where(z>1100,sound_speed(z,"photons"),sound_speed(z,type_))
    return c*np.sqrt(np.pi/(G*density(z)))

def jeans_mass(z):
    G = constants.G
    rho = density(z)
    c = np.where(z>1100,sound_speed(z,"photons"),sound_speed(z,"gas"))
    return (np.pi)**(2./5)/6.*c**3/(G**(3./2)*rho**(1/2.))


if __name__ == '__main__':
    #a)
    """
    a = np.linspace(1e-4,1,1000)
    plt.loglog(a,T_photon(a),label="Photons")
    plt.loglog(a,T_gas(a),label="Gas")
    plt.title("Temperature of gas and photons",fontsize=20)
    plt.xlabel("a",fontsize=15)
    plt.ylabel("T [K]",fontsize=15)
    plt.legend(loc="best",fontsize=15)
    plt.show()
    """
    #b)
    z = np.linspace(2000,0,10000)
    jeans_length_gas = jeans_length(z,"gas")
    jeans_mass = jeans_mass(z)
    plt.loglog(z,jeans_length_gas,label="Gas")
    plt.legend()
    plt.xlim(z[0],z[-1])
    plt.show()
    plt.loglog(z,jeans_mass)

    plt.xlim(z[0],z[-1])
    plt.show()
