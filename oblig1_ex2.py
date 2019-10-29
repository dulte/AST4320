import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import constants

try:
    import seaborn as sns
    sns.set_style("darkgrid")
except:
    print("Everything Looks Nicer With Seaborn!")


rc('font',**{'family':'serif'})


def H(omega_m,omega_l,H0,a):
    return H0*np.sqrt(omega_m*a**(-3)+omega_l)

def a_dot(omega_m,omega_l,H0,a):
    return H0*np.sqrt(omega_m*a**(-1)+omega_l*a**2)


def dxdt(omega_m,omega_l,a,H0,delta,x):

    H_ = H(omega_m,omega_l,H0,a)
    return delta*3/2.*H0**2*(omega_m*a**(-3)+omega_l) - 2*H_*x


def solve(omega_m,omega_l,delta0,a0,steps,H0):

    a_array = np.linspace(a0,1,steps)
    delta_array = np.zeros_like(a_array)
    delta_array[0] = delta0
    x_array = np.zeros_like(a_array)
    x_array[0] = a_dot(omega_m,omega_l,H0,a_array[0])
    da = abs(a_array[1]-a_array[0])


    for i in range(steps-1):
        a_d = a_dot(omega_m,omega_l,H0,a_array[i])
        dx = dxdt(omega_m,omega_l,a_array[i],H0,delta_array[i],x_array[i])*da*1./a_d
        x_array[i+1] = x_array[i] + dx
        ddelta = x_array[i+1]*da*1./a_d
        delta_array[i+1] = delta_array[i] + ddelta
    

    return delta_array,a_array, x_array

def getF(delta,a):
    return (np.log(delta[1:]/delta[:-1]))/(np.log(a[1:]/a[:-1]))
    


if __name__ == '__main__':

    H0 = 2e-18
    delta0 = 1e-3
    a0 = 1e-3

    steps = int(1e5)


    omega_ms = [1,0.3,0.8]
    fs = []
    a_sols = []
    for omega_m in omega_ms:
        omega_l = 1- omega_m
        print(omega_m,omega_l)
        delta_sol,a_sol,delta_dot_sol = solve(omega_m,omega_l,delta0,a0,steps,H0)

        a_sols.append(a_sol)
        fs.append(getF(delta_sol,a_sol))
        plt.loglog(a_sol,delta_sol,label=r"$\Omega_m=$%.1f,$\Omega_{\Lambda}=$%.1f"%(omega_m,omega_l))
        print("Done with omega_m=%s"%omega_m)

    plt.title(r"Density Contrast $\delta$ for Three Different Universes",fontsize=20)
    plt.xlabel("a",fontsize=15)
    plt.ylabel(r"$\delta$",fontsize=15)
    plt.legend(loc="best",fontsize=15)
    plt.show()


    for i,f in enumerate(fs):

        omega_m = omega_ms[i]
        omega_l = 1-omega_m
        z = 1./a_sols[i][1:] -1
        plt.loglog(z,f,label=r"$\Omega_m=$%.1f,$\Omega_{\Lambda}=$%.1f"%(omega_m,omega_l))

    plt.title(r"Growth Factor $f$ for Three Different Universes",fontsize=20)
    plt.xlabel("z",fontsize=15)
    plt.ylabel(r"f",fontsize=15)
    plt.legend(loc="best",fontsize=15)
    plt.ylim(0,1.1)
    plt.xlim(z[0],z[-1])
    plt.show()
