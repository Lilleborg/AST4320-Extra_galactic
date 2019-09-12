import numpy as np 
import matplotlib.pyplot as plt
import astropy.units as unit
import matplotlib as mpl

# Plotting style
plt.style.use('bmh')

def Hubbleparam(a,Omega_m,Omega_lam):
    res = H_0*np.sqrt((Omega_m*a**(-3)+Omega_lam))
    return res

def dx_da(x,delta,a,cosmo):
    # cosmo = [0.8,0.2] ex, array containing fraction of each component
    res = 3./2*Hubbleparam(a,*cosmo)*delta/a-2*x/a
    return res

# Constants and parameters
H_0 = 70# * unit.km / unit.s / unit.Mpc
cosmos = [[1.0,0.0],[0.3,0.7],[0.8,0.2]]

# Integration variable, scale factor
steps = int(1e4)
a_0 = 1
a_start = 10**(-3)
a = np.linspace(a_start,a_0,steps)
z = 1/a - 1
da = np.abs(a[1]-a[0])

delta = np.zeros_like(a)
delta[0] = a_start
delta_dot = np.zeros_like(delta)

fig, ax = plt.subplots(2,1)

for cosmo in cosmos:
    # Boundary for derivative dependent on cosmology
    delta_dot[0] = Hubbleparam(a[0],*cosmo)*a[0]

    for i in range(len(a)-1):
        delta_dot[i+1] = delta_dot[i] + dx_da(delta_dot[i],delta[i],a[i],cosmo)*da
        delta[i+1] = delta[i] + delta_dot[i+1]/Hubbleparam(a[i],*cosmo)/a[i]*da
        
    f = delta_dot/delta/Hubbleparam(a,*cosmo)

    ax[0].loglog(a,delta)
    ax[1].loglog(z,f,label=r'$(\Omega_m,\Omega_\Lambda) = (%.1f,%.1f)$'%(cosmo[0],cosmo[1]))

# h, l = ax[1].get_legend_handles_labels()
# ax[0].legend(h,l)
ax[1].legend()

fig.suptitle('Time evolution of perturbations in different cosmologies')
ax[0].set_title('Overdensity vs scale factor')
ax[0].set_xlabel(r'$a$')
ax[0].set_ylabel(r'$\delta$')

ax[1].set_title('Growth factor vs red shift')
ax[1].set_xlabel(r'$z$')
ax[1].set_ylabel(r'$f$')
ax[1].set_xlim(z[0]+450,z[-2])
# ax[1].use_sticky_edges = False
# ax[1].margins(x=20,y=0.05,tight=None)

plt.show()
fig.savefig('exe2.pdf')