import numpy as np
import matplotlib.pyplot as plt

plt.style.use('bmh')

R = 1
k = np.linspace(-20,20,10000)

def w(k,R):
    res = np.where(k==0,2*R,2*np.sin(k*R)/k)
    fwhm = FWHM(res,k)
    return res,fwhm

def FWHM(w,k):
    half_max = w.max()/2
    id_half_max = np.argmin(np.abs(w-half_max)) # first id at half max
                                        # this gives k value on negative side
    k_value = np.abs(k[id_half_max])
    return [[-k_value,k_value],half_max]

w,fwhm = w(k,R)
width = np.abs(fwhm[0][1]-fwhm[0][0])
print('The full width at half max: FWHM = {:.2f}'.format(width))

fig,ax = plt.subplots()
ax.plot(k,w)
#ax.annotate(s='FWHM',xy=(fwhm[0][0],fwhm[1]),xytext=(fwhm[0][1],fwhm[1]),arrowprops=dict(arrowstyle='<->'))
ax.arrow(fwhm[0][0],fwhm[1],width,0,width=0.01)# not so fancy, but does the trick
ax.set_xlabel('k')
ax.set_ylabel(r'$\tilde{W}(k)$')
ax.set_title('Fourier Transformed Top-hat Smoothing Function in 1D')
fig.savefig('window_func.pdf')
plt.show()