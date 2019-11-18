import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcdefaults()
plt.style.use('bmh')

sigma_T = 6.65e-25  # cm^2
c = 3e10            # cm/s
H_0 = 2.2e-18       # s^-1
n_e = 1.9e-7        # cm^-3
Omega_L = 0.692
Omega_m = 0.308
zs = np.linspace(0,10,1000)
taus = np.zeros_like(zs)

def integrand(z):
    prefactor = c*sigma_T*n_e/H_0
    return prefactor*(1+z)**2/(np.sqrt(Omega_L+Omega_m*(1+z)**3))

for i,z in enumerate(zs):
    taus[i] = integrate.quad(integrand,0,z)[0]

intresting_zs = [0,6,10]
for z in intresting_zs:
    print('Optical depth at z = {}, tau = {:.4f}'.format(z,taus[np.argmin(np.abs(zs-z))]))

fig, ax = plt.subplots()
ax.plot(zs,taus)
ax.set_title('Optical depth of ionized IGM')
ax.set_ylabel(r'$\tau_e(z)$')
ax.set_xlabel(r'$z$')
fig.savefig('optical_depth_exe2.pdf')
plt.show()