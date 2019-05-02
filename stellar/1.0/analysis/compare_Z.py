

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









# Zs = [0.001, 0.01]
Zs = [0.00001, 0.0001, 0.001, 0.01]

age = 6.0

SPS = 'BPASSv2.1.binary/ModSalpeter_100'

for ia,Z in enumerate(Zs):

    c = cm.plasma(float(ia)/float(len(Zs)+1))

    lam = np.load('../inputs/outputs/'+SPS.split('/')[0]+'/lam.npy')
    
    Lnu = np.load('../inputs/outputs/'+SPS+'/'+str(Z)+'_'+str(age)+'.npy')

    # if SPS == 'P2/ModSalpeter_100': Lnu /= 1E6

    s = [(lam<10000)&(lam>250)]

    ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = Z, c = c, lw = 1)
    










ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/Z.pdf')