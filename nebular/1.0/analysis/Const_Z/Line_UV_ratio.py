

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle
from scipy import integrate

plt.style.use('simple')



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15  
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI4861','OIII4959','OIII5007']

labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'\regular{[CII]1907,1909}'
labels['OII3726-OII3729'] = r'\regular{[OII]3726,3729}'
labels['NeIII3869'] = r'\regular{[NeIII]3869}'
labels['NeIII3967'] = r'\regular{[NeIII]3967}'
labels['HI4340'] = r'\regular{H\gamma}'
labels['HI4861'] = r'\regular{H\beta}'
labels['OIII5007'] = r'\regular{[OIII]5007}'
labels['OIII4959'] = r'\regular{[OIII]4959}'
labels['HI6563'] = r'\regular{H\alpha}'



continuum_grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/nebular.p','rb'), encoding='latin1')

line_grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')



s = np.fabs(continuum_grid['lam']-1500.)<200.
lam = continuum_grid['lam'][s]

log10tmaxs = np.arange(7.,9.1,0.1)

for iZ, Z, ls in zip(range(len(continuum_grid['Z'])), continuum_grid['Z'], ['-','-.','--',':']*3):

    c = cm.viridis(float(iZ)/float(len(continuum_grid['Z'])+1))

  
    LUV = []
    Lline = []
   
    for log10tmax in log10tmaxs:
    
        SED = np.zeros(len(continuum_grid['lam']))
        
        line_luminosity = 0.0

        for ia, log10age in enumerate(np.arange(6.,log10tmax+0.1,0.1)):

            # --- determine SF activity in bin

            SF = 10**log10age - 10**(log10age-0.1) 
            if log10age==6.0: SF = 10**6.0     

            SED += SF * (continuum_grid['stellar'][ia,iZ] + continuum_grid['nebular'][ia,iZ])
    
            for l in line:
            
                l = l.encode('UTF-8')
            
                line_luminosity += SF*10**line_grid[l]['luminosity'][ia, iZ]

        LUV.append(integrate.trapz((1./lam) * SED[s], x = lam) / integrate.trapz((1./lam), x = lam))
        

        Lline.append(line_luminosity)

    LUV = np.log10(np.array(LUV)) + np.log10(3E8/1500E-10)
    Lline = np.log10(np.array(Lline))

    R = Lline - LUV

    ax.plot(log10tmaxs - 6., R, label = r'$\rm Z='+str(Z)+r'$', c=c, ls=ls, zorder = 1, lw = 1)
    




        
# ax.legend(loc = 'center left', fontsize=7, bbox_to_anchor = (1.0, 0.5))

ax.legend(loc = 'lower left', fontsize=6, handlelength=2.5, labelspacing=0.2)


ax.set_xlim([1.,3.])
# ax.set_ylim([-3.5, 3.6])


ax.set_xlabel(r'$\rm\log_{10}(t_{SF}/Myr)$')
ax.set_ylabel(r'$\rm\log_{10}(L_{[OIII]+H\beta}/L^{FUV})$')


figname = 'figures/Line_FUV.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



