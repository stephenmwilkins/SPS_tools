

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
ages = [6.0, 7.0, 8.0, 9.0, 10.0]


# SPSs = ['BPASSv2.1.binary/ModSalpeter_100', 'BPASSv2.1.binary/ModSalpeter_300']

SPSs = ['BPASSv2.1.binary/ModSalpeter_100', 'P2/ModSalpeter_100']

SPSs = ['BPASSv2.1.binary/ModSalpeter_100', 'BPASSv2.binary/ModSalpeter_100']

for ia,age in enumerate(ages):

    c = cm.plasma(float(ia)/float(len(ages)+1))

    for SPS, ls in zip(SPSs, ['-','-.',':']):
    
        lam = np.load('../inputs/outputs/'+SPS.split('/')[0]+'/lam.npy')
        
        Lnu = np.load('../inputs/outputs/'+SPS+'/'+str(Z)+'_'+str(age)+'.npy')

        # if SPS == 'P2/ModSalpeter_100': Lnu /= 1E6

        s = [(lam<10000)&(lam>100)]

        ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = SPS, c = c, ls = ls)
        










ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/SPS.pdf')