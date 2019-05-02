import numpy as np
import pickle
import os
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.style.use('simple')



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15  
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


# --- Build stellar SED assuming continuous star formation activity.


# SPS = 'P2'
# IMF = 'ModSalpeter_100'

SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/nebular.p','rb'), encoding='latin1')

print(grid.keys())
print(grid['stellar'].shape)
print(grid['stellar'][0,4][300:320])


s = np.fabs(grid['lam']-1500.)<200.
lam = grid['lam'][s]


log10tmaxs = np.arange(7.,9.1,0.1)

for iZ, Z, ls in zip(range(len(grid['Z'])), grid['Z'], ['-','-.','--',':']*3):

    c = cm.viridis(float(iZ)/float(len(grid['Z'])+1))

    LUV = []

    for log10tmax in log10tmaxs:
    
        SED = np.zeros(len(grid['lam']))

        for ia, log10age in enumerate(np.arange(6.,log10tmax+0.1,0.1)):

            # --- determine SF activity in bin

            SF = 10**log10age - 10**(log10age-0.1) 
            if log10age==6.0: SF = 10**6.0     

            SED += SF * (grid['stellar'][ia,iZ] + grid['nebular'][ia,iZ])
            # SED += SF * (grid['stellar'][ia,iZ])

        luv =  integrate.trapz((1./lam) * SED[s], x = lam) / integrate.trapz((1./lam), x = lam)

        # LUV.append(np.interp(1500., grid['lam'], SED))
        LUV.append(luv)
        


    ax.plot(log10tmaxs - 6., np.log10(np.array(LUV)), label = r'$\rm Z='+str(Z)+r'$', c=c, ls=ls)


   
   
        
# ax.legend(loc = 'center left', fontsize=7, bbox_to_anchor = (1.0, 0.5))

ax.legend(loc = 'upper left', fontsize=6, handlelength=2.5, labelspacing=0.2)


ax.set_xlim([1.,3.])
ax.set_ylim([27.8, 28.3])


ax.set_xlabel(r'$\rm\log_{10}(t_{SF}/Myr)$')
ax.set_ylabel(r'$\rm\log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')



C = 43.35

B = C - np.log10(3E8/1500E-10)

print(B) 

ax.axhline(B, c='k', alpha=0.1, lw=3, zorder = -1)
ax.text(2., B+0.005, r'Kennicutt & Evans (2010)', fontsize = 6, color='k', zorder = -1, alpha = 0.5)


figname = 'figures/SFR_UV_age.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot