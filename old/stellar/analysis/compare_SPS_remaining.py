

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate




plt.style.use('simple')

fig = plt.figure( figsize=(6,4) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))









Z = 0.004
Z = 0.008
# Z = 0.02


SPSs = ['BPASSv2.1.single/ModSalpeter_100','BPASSv2.1.binary/ModSalpeter_100', 'BPASSv2.1.binary/ModSalpeter_300', 'P2/ModSalpeter_100']

for SPSIMF, ls in zip(SPSs, ['-','--','-.',':']):
 
    mass = np.load('../SSP/outputs/'+SPSIMF+'/mass_'+str(Z)+'.npy')
 
    
 
    ages = np.indices(mass.shape)[0]*0.1 + 6.
 
    print SPSIMF, len(mass), len(ages)

    ax.plot(ages, mass, label = SPSIMF, c='k', ls = ls)
        










ax.legend()

ax.set_ylabel(r'$f_{\rm rem}$')
ax.set_xlabel(r'$\log_{10}(age/yr)$')
    
fig.savefig('figures/SPS_remaining.pdf')