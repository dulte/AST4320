import numpy as np
import matplotlib.pyplot as plt
from scipy import constants,stats
from scipy.integrate import ode

def H(omega_m,omega_l,t,H0):
    if omega_m == 1:
        return 2/(3*t)
    else:
        return H0*np.sqrt(omega_l)*np.cosh(3./2*np.sqrt(omega_l)*H0*t)\
                            *(np.sinh(3./2*np.sqrt(omega_l)*H0*t))**(1/3)

def t(omega_m,omega_l,a,H0):
    if omega_m == 1:
        return 2./3*1/H0*a**(3./2)
    else:
        return 2./3*1/(np.sqrt(omega_l)*H0)*np.arcsinh(((omega_m/omega_l)**(-1./3)*a)**(3./2))

def a_dot(omega_m,omega_l,t,H0):
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
    H_ = H(omega_m,omega_l,t,H0)
    rho0 = 3*H0**2/(8*np.pi*G)
    return delta*4*np.pi*rho0*a**(-3) - 2*H_*x

def ddeltadt(x):
    return x

def delta_step(t,delta,args):
    x = args[0]
    return ddeltadt(x)
def x_step(t,delta,args):
    omega_m,omega_l,a,H0,x = args[0],args[1], args[2],args[3],args[4]
    return dxdt(omega_m,omega_l,a,t,H0,delta,x)


def solve(omega_m,omega_l,delta0,a0,steps,H0):
    t_array = np.linspace(t(omega_m,omega_l,a0,H0),t(omega_m,omega_l,1,H0),steps)
    #print(t_array)
    #exit()
    a_array = a(omega_m,omega_l,t_array,H0)
    delta_array = np.zeros_like(a_array)
    delta_array[0] = delta0
    x_array = np.zeros_like(a_array)
    x_array[0] = a_dot(omega_m,omega_l,t_array[0],H0)

    dt = t_array[1]-t_array[0]


    r_delta = ode(delta_step).set_integrator('lsoda', method='bdf')
    r_delta.set_initial_value(delta_array[0],t_array[0]).set_f_params([x_array[0]])

    r_x = ode(x_step).set_integrator('lsoda', method='bdf')
    r_x.set_initial_value(x_array[0],t_array[0]).set_f_params([omega_m,omega_l,a_array[0],H0,x_array[0]])
    """
    for i in range(steps-1):

        dx = dxdt(omega_m,omega_l,a_array[i],t_array[i],H0,delta_array[i],x_array[i])*dt
        x_array[i+1] = x_array[i] + dx
        ddelta = x_array[i+1]*dt
        delta_array[i+1] = delta_array[i] + ddelta
        if delta_array[i+1] != delta_array[i+1]:
            break
    """
    i = 0
    while r_delta.successful() and r_x.successful() and r_delta.t <= t_array[-1]:
        r_x.set_f_params([omega_m,omega_l,a_array[i],H0,x_array[i]])

        x_array[i+1] = r_x.integrate(r_x.t + dt)[0]

        r_delta.set_f_params([x_array[i+1]])
        delta_array[i+1] = r_delta.integrate(r_delta.t + dt)[0]


    return delta_array,a_array

def find_power(delta,a):
    pass

if __name__ == '__main__':
    omega_m = 1.
    omega_l = 0.
    H0 = 2e-18
    delta0 = 1e-3
    a0 = 1e-3

    steps = int(1e6)

    delta,a = solve(omega_m,omega_l,delta0,a0,steps,H0)
    print(delta)
    plt.loglog(a,delta)
    plt.show()
