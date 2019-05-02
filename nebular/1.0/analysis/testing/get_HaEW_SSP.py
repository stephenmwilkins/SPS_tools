



import numpy as np
import cPickle as pickle
from scipy import interpolate
from scipy import integrate
import read

import matplotlib.pyplot as plt
from matplotlib import cm

c = 3E8


SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_100'

# SPS = 'P2'
# IMF = 'ModSalpeter_100'




# --- define plot 

plt.style.use('simple')

fig = plt.figure( figsize=(5,3) )

left  = 0.1 
bottom = 0.1   
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))








ages = np.arange(6.0,9.1,0.1)

Zs = [0.00001,0.0001,0.001,0.01,0.04]
# Zs = [0.001]

for include_nebular, ls in zip([True, False],['-','-.']):
  

    for iZ,Z in enumerate(Zs):

        c = cm.viridis(float(iZ)/float(len(Zs)+1))

        EWs = []
    
        print Z

        for log10age in ages:
    
            print log10age

            model = str(Z)+'_'+str(log10age)

            lines = ['HI6563'] #,'HI4861','CIII1909','CIII1907'

            lines = np.array(lines)

            lam_lines, emergent  = read.lines(SPS, IMF, model, lines = lines)

            lam_continuum, nebular, nebular_wlines = read.nebular_continuum(SPS, IMF, model)

            lam_continuum, stellar = read.stellar_continuum(SPS, IMF, model)

            if include_nebular:

                total = stellar  + nebular
                
            else:
            
                total = stellar

            lam = lam_lines[0]

            EWs.append(lam*(10**emergent[0])/np.interp(lam, lam_continuum[::-1], total[::-1]))


        ax.plot(ages, np.log10(np.array(EWs)), c=c, label='Z='+str(Z), lw = 0.5, ls = ls) 
    

        
        
        
        
        
        
        
ax.legend(fontsize=6, labelspacing = 0.0)



ax.set_ylabel(r'$\rm\log_{10}(EW_{H\alpha}/\AA)$')
ax.set_xlabel(r'$\rm \log_{10}(age/yr)$')
    
fig.savefig('figures/HaEW.pdf')



