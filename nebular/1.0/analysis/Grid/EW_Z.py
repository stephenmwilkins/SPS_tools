

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import readGrid
import copy





SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'




line = 'HI1216'
# line = 'CIII1907,CIII1909'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = True)

# Zs = grid.nebular_line_grid['Z']
# log10Zs = grid.nebular_line_grid['log10Z']

Zs = np.array([0.00001, 0.0001, 0.001, 0.004, 0.01, 0.04])
log10Zs = np.log10(Zs)

log10ages = grid.stellar_grid['log10age']

# ---- define default parameters

fesc = 0.0

log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}



# ---- define plot


plt.style.use('simple')

fig = plt.figure( figsize=(3,3) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))




for iZ, log10Z in enumerate(log10Zs):

    p = copy.copy(default_p)

    p['log10Z'] = log10Z

    EW = []
    
    c = cm.plasma(float(iZ)/(len(log10Zs)-1))
    
    for log10age in log10ages:
    
        p['log10age'] = log10age
    
        grid.generate_nebular_continuum(p)

        grid.generate_stellar(p)
    
        EW.append(np.log10(grid.EW(p, line)))

    precision = int(np.ceil(abs(np.log10(abs(Zs[iZ]))))) if p != 0 else 1
    print log10Z, precision
    ax.plot(log10ages, EW, label = '$'+'{:.{}f}'.format(Zs[iZ], precision)+'$', c = c, alpha = 0.8, lw =0.5)


# model_label = ' \quad '.join([SPS_name, 'IMF='+IMF, 'Z='+str(p['Z']), r'\log_{10}(n/cm^{-3})='+str(int(np.round(p['log10n_H'],0))), '}$\n ${\\rm',r'\log_{10}(U_S)='+str(int(np.round(p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(p['CO']), r'\xi_{d}='+str(p['d2m'])])

model_label = ' \quad '.join([ r'\log_{10}(n/cm^{-3})='+str(int(np.round(p['log10n_H'],0))), r'\log_{10}(U_S)='+str(int(np.round(p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(p['CO']), r'\xi_{d}='+str(p['d2m'])])

print model_label


plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 5, alpha = 0.6, ha = 'center', va = 'bottom')


ax.legend(prop={'size':6}, labelspacing = 0.0, title = r'$\rm\log_{10}(Z)$')
ax.legend(prop={'size':6}, labelspacing = 0.0)


ax.set_ylabel(r'$\rm\log_{10}(EW/\AA)$')
ax.set_xlabel(r'$\rm\log_{10}(t_{SF}/yr)$')
    
fig.savefig('figures/EW_Z_'+line+'.pdf')



