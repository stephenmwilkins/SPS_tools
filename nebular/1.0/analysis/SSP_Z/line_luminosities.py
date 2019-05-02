

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.1
bottom = 0.1  
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))



ax.axhline(0.0,c='k',alpha=0.1,lw=2)

ax.fill_between([-3.0,-2.0], [-7,-7], [2.,2.], color='k',alpha=0.05)



SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI6563']
# line = 'CIII1907,CIII1909'

lines = [['HI4861'],['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563']]


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




grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')

print(grid.keys())


for il, line, ls in zip(range(len(lines)), lines, ['-','-.','--',':']*3):

    c = cm.magma(float(il-1)/float(len(lines)))

    L = []
    
    label = '-'.join(line)

    for iZ, Z in enumerate(grid['Z']):

        line_luminosity = 0.0

        for l in line:
        
            l = l.encode('UTF-8')
        
            line_luminosity += 10**grid[l.encode('UTF-8')]['luminosity'][0, iZ]

        L.append(line_luminosity)
        
    L = np.log10(np.array(L))
    
    if il == 0: 
        L_Hbeta = L    
    if line[0] != 'HI4861':
        ax.plot(grid['log10Z'], L - L_Hbeta, label = r'$\rm '+labels[label]+r'$', c=c, ls=ls, lw=1)


        
# ax.legend(loc = 'center left', fontsize=6, bbox_to_anchor = (1.0, 0.5))
ax.legend(loc = 'lower left', fontsize=6, handlelength=2.5, labelspacing=0.2)


ax.set_xlim([-5.0, -1.4])
ax.set_ylim([-6.5,1.9])

ax.set_xlabel(r'$\rm\log_{10}(Z)$')
ax.set_ylabel(r'$\rm\log_{10}[(L/L_{H\beta})_{1\ Myr}]$')


figname = 'figures/line_luminosities.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



