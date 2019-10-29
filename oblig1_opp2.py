import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import constants

rc('font',**{'family':'serif'})


def H(omega_m,omega_l,t,H0,a):

    return H0*np.sqrt(omega_m*a**(-3)+omega_l)

    if omega_m == 1:
        return 2/(3*t)
    else:
        return H0*np.sqrt(omega_l)*np.cosh(3./2*np.sqrt(omega_l)*H0*t)\
                            *(np.sinh(3./2*np.sqrt(omega_l)*H0*t))**(-1)

def t(omega_m,omega_l,a,H0):
    if omega_m == 1:
        return 2./3*1/H0*a**(3./2)
    else:
        return 2./3*1/(np.sqrt(omega_l)*H0)*np.arcsinh(((omega_m/omega_l)**(-1./3)*a)**(3./2))

def a_dot(omega_m,omega_l,t,H0,a):
    return H0*np.sqrt(omega_m*a**(-1)+omega_l*a**(2))
    if omega_m == 1:
        return (3./2*H0)**(2./3)*t**(-1./3)*2./3
    else:
        return (omega_m/omega_l)**(1./3)*(np.sinh(3./2*np.sqrt(omega_l)*H0*t))**(-1/3)*\
                np.cosh(3./2*np.sqrt(omega_l)*H0*t)*np.sqrt(omega_l)*H0

def a(omega_m,omega_l,t,H0):
    if omega_m == 1:
        return (3./2*H0*t)**(2./3)
    else:
        return (omega_m/omega_l)**(1./3)*(np.sinh(3./2*np.sqrt(omega_l)*H0*t))**(2/3)


def dxdt(omega_m,omega_l,a,t,H0,delta,x):
    G = constants.G
    H_ = H(omega_m,omega_l,t,H0,a)
    return delta*3/2.*H0**2*a**(-3)- 2*H_*x
    rho0 = 3*H0**2/(8*np.pi*G)
    return delta*4*np.pi*rho0*a**(-3)*G - 2*H_*x

def solve(omega_m,omega_l,delta0,a0,steps,H0):
    t_array = np.linspace(t(omega_m,omega_l,a0,H0),t(omega_m,omega_l,1,H0),steps)
    a_array = a(omega_m,omega_l,t_array,H0)
    delta_array = np.zeros_like(a_array)
    delta_array[0] = delta0
    x_array = np.zeros_like(a_array)
    x_array[0] = a_dot(omega_m,omega_l,t_array[0],H0,a_array[0])
    dt = t_array[1]-t_array[0]


    for i in range(steps-1):

        dx = dxdt(omega_m,omega_l,a_array[i],t_array[i],H0,delta_array[i],x_array[i])*dt
        x_array[i+1] = x_array[i] + dx
        ddelta = x_array[i+1]*dt
        delta_array[i+1] = delta_array[i] + ddelta
        if delta_array[i+1] != delta_array[i+1]:
            break

    return delta_array,a_array


if __name__ == '__main__':

    H0 = 2e-18
    delta0 = 1e-3
    a0 = 1e-3

    steps = int(1e5)


    omega_ms = [1,0.3,0.8]

    for omega_m in omega_ms:
        omega_l = 1- omega_m
        print(omega_m,omega_l)
        delta_sol,a_sol = solve(omega_m,omega_l,delta0,a0,steps,H0)

        plt.loglog(a_sol,delta_sol,label=r"$\Omega_m=$%.1f,$\Omega_{\Lambda}=$%.1f"%(omega_m,omega_l))
        print("Done with omega_m=%s"%omega_m)

    plt.title(r"Density Contrast $\delta$ for Three Different Universes",fontsize=20)
    plt.xlabel("a",fontsize=15)
    plt.ylabel(r"$\delta$",fontsize=15)
    plt.legend(loc="best")
    plt.show()
