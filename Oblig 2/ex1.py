import numpy as np
import matplotlib.pyplot as plt


def W(k,R):
    return np.where(k == 0,2*R,2*np.sin(k*R)/k)

def find_FWHM(k,W):
    W_max = np.max(W)
    HM = W_max/2.
    FWHM = k[np.argmin(np.abs(W-HM))]
    return [[FWHM,-FWHM],HM]


k = np.linspace(-8*np.pi,8*np.pi,1000)
R = 1

W_array = W(k,R)

FWHM = find_FWHM(k,W_array)
length = np.max(FWHM[0]) - np.min(FWHM[0])



print("FWHM = %g"%(length))
plt.plot(k,W_array,label=r"$\tilde{W}(k)$")
plt.annotate(s='', xy=(np.max(FWHM[0]),FWHM[1]), xytext=(np.min(FWHM[0]),FWHM[1]), arrowprops=dict(arrowstyle='<->'),label="FWHM")
#plt.legend()
plt.xlabel(r"$k$",fontsize=25)
plt.ylabel(r"$\tilde{W}(k)$",fontsize=25)
plt.title("Fourier Transform of 1D Top-hat Smoothing Function",fontsize=25)
plt.show()
