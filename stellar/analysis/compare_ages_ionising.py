

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

ages = np.arange(6.4,6.9,0.1)

SPS = 'BPASSv2.1.binary/ModSalpeter_100'

for ia,age in enumerate(ages):

    c = cm.plasma(float(ia)/float(len(ages)+1))

    lam = np.load('../SSP/outputs/'+SPS.split('/')[0]+'/lam.npy')
    
    Lnu = np.load('../SSP/outputs/'+SPS+'/'+str(Z)+'_'+str(age)+'.npy')

    s = [(lam<10000)&(lam>250)]

    ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = age, c = c, lw = 1)
    










ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/ages.pdf')