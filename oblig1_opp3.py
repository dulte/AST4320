import numpy as np
import matplotlib.pyplot as plt
try:
    import seaborn as sns
    sns.set_style("darkgrid")
except:
    print("Everything Looks Nicer With Seaborn!")


def T_photon(a):
    return 4*1./a
def T_gas(a):
    return 4e-3*1./a**(2)

if __name__ == '__main__':
    a = np.linspace(1e-4,1,1000)
    plt.loglog(a,T_photon(a),label="Photons")
    plt.loglog(a,T_gas(a),label="Gas")
    plt.title("Temperature of gas and photons",fontsize=20)
    plt.xlabel("a",fontsize=15)
    plt.ylabel("T [K]",fontsize=15)
    plt.legend(loc="best",fontsize=15)
    plt.show()
