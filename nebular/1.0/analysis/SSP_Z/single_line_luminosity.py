

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')



fig = plt.figure(figsize=(3.5,3.5))

left  = 0.15  
bottom = 0.15   
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))




SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI6563']
# line = 'CIII1907,CIII1909'

labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'C\regular_{III}'
labels['OII3726-OII3729'] = r'O\regular_{II}'
labels['NeIII3869'] = r'Ne\regular_{III}'
labels['HI4861'] = r'H\beta'
labels['OIII4363-OIII5007'] = r'O\regular_{III}'
labels['HI6563'] = r'H\alpha'


grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')

print(grid.keys())


for iZ, Z, ls in zip(range(len(grid['Z'])), grid['Z'], ['-','-.','--',':']*3):

    print(Z)

    color = cm.viridis(float(iZ)/float(len(grid['Z'])+1))

    L = []

    for ia,log10age in enumerate(grid['log10age']):

        lum = 0.0

        for l in line:

            l = l.encode('UTF-8')

            lum += grid[l]['luminosity'][ia, iZ]
        
        L.append(lum)
        
    if iZ==0: 
        Z_label = '0.00001' 
    else:
        Z_label = str(Z)
    
    ax.plot(grid['log10age'], np.array(L), label = r'$\rm Z='+Z_label+'$', c=color, ls=ls, lw=1)
        

ax.legend(fontsize=6)
ax.set_xlim([6.,9.])

label = '-'.join(line)

ax.set_xlabel(r'$\rm\log_{10}(age/yr)$')
ax.set_ylabel(r'$\rm\log_{10}(L_{'+labels[label]+'}/erg\,s^{-1}\,M_{\odot}^{-1})$')

fig.savefig('figures/line_luminosities_SSP_'+str(label)+'.pdf')

