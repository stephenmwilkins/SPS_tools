



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




lines = ['HI6563'] #,'HI4861','CIII1909','CIII1907'
lines = np.array(lines)



log10ages = np.arange(6.0,9.1,0.1)

Zs = [0.00001,0.0001,0.001,0.01,0.04]
# Zs = [0.001]



for iZ,Z in enumerate(Zs):

    c = cm.viridis(float(iZ)/float(len(Zs)+1))

    EWs_stellar = []
    EWs_total = []

    totSF = 0.0

    for i,log10age in enumerate(log10ages):

        SF = 10**log10age - totSF

        totSF += SF

        model = str(Z)+'_'+str(log10age)

        lam_lines, emergent  = read.lines(SPS, IMF, model, lines = lines)
        lam_continuum, nebular, nebular_wlines = read.nebular_continuum(SPS, IMF, model)
        lam_continuum, stellar = read.stellar_continuum(SPS, IMF, model)

        # --- include nebular continuum emission

        total = stellar  + nebular

        # --- add up

        if i==0:
            line_total = SF * 10**emergent 
            total_total = SF * total
            stellar_total = SF * stellar
        else:
            line_total += SF * 10**emergent
            total_total += SF * total  
            stellar_total += SF * stellar
    


        lam = lam_lines[0]
        EW_total = lam*(line_total[0])/np.interp(lam, lam_continuum[::-1], total_total[::-1])
        EW_stellar = lam*(line_total[0])/np.interp(lam, lam_continuum[::-1], stellar_total[::-1])
        
        print log10age, EW_total, EW_stellar

        EWs_total.append(EW_total)
        EWs_stellar.append(EW_stellar)


    ax.plot(log10ages[log10ages>6.9], np.log10(np.array(EWs_total)[log10ages>6.9]), c=c, label='Z='+str(Z), lw = 0.5, ls = '-') 
    ax.plot(log10ages[log10ages>6.9], np.log10(np.array(EWs_stellar)[log10ages>6.9]), c=c, lw = 0.5, ls = '-.') 


    
        
        
        
        
        
        
ax.legend(fontsize=6, labelspacing = 0.0)



ax.set_ylabel(r'$\rm\log_{10}(EW_{H\alpha}/\AA)$')
ax.set_xlabel(r'$\rm \log_{10}(age/yr)$')
    
fig.savefig('figures/HaEW_const.pdf')



