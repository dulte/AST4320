import numpy as np
import matplotlib.pyplot as plt


def T_photon(a):
    return 1./a
def T_gas(a):
    return 1./a**(2)

if __name__ == '__main__':
    a = np.linspace(1e-4,1,1000)
    plt.loglog(a,T_photon(a),label="Photons")
    plt.loglog(a,T_gas(a),label="Gas")
    plt.title("Temperature of gas and photons")
    plt.xlabel("a")
    plt.ylabel("T [K]")
    plt.legend()
    plt.show()
