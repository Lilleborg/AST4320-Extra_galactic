import numpy as np 
import matplotlib.pyplot as plt
import astropy.units as unit
import matplotlib as mpl

# Plotting style
plt.style.use('bmh')

a = np.logspace(-4,0,1e4)

T_gamma = 2.725/a 
T_gas = 2.725/1091/a**2

fig,ax = plt.subplots()

ax.semilogy(a,T_gamma,label=r'$T_\gamma^{bb}$')
ax.semilogy(a,T_gas,label=r'$T_{gas}$')

fig.suptitle('Time evolution of gas and radiation temperatures')
ax.legend()
ax.set_xlabel(r'$a$')
ax.set_ylabel(r'$T$')
fig.savefig('exe3.pdf')
plt.show()