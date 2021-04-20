

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




SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = ['HI6563']
# line = 'CIII1907,CIII1909'

label = '-'.join(line)

labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'C\regular_{III}'
labels['OII3726-OII3729'] = r'O\regular_{II}'
labels['NeIII3869'] = r'Ne\regular_{III}'
labels['HI4861'] = r'H\beta'
labels['OIII4363-OIII5007'] = r'O\regular_{III}'
labels['HI6563'] = r'H\alpha'


grid = pickle.load(open('/Users/stephenwilkins/research/data/SPS/nebular/2.0/Z_refQ/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')

print(grid.keys())




for iZ, Z, ls in zip(range(len(grid['Z'])), grid['Z'], ['-','-.','--',':']*3):

    color = cm.viridis(float(iZ)/float(len(grid['Z'])+1))

    EWs = []

    for ia,log10age in enumerate(grid['log10age']):

        line_luminosity = 0.0
        stellar_continuum = 0.0
        nebular_continuum = 0.0
        lam = 0.0

        for l in line:

            # l = l.encode('UTF-8')

            line_luminosity += 10**grid[l]['luminosity'][ia, iZ]
            stellar_continuum += grid[l]['stellar_continuum'][ia, iZ]
            nebular_continuum += grid[l]['nebular_continuum'][ia, iZ]
            lam += grid[l]['lam']

        total_continuum = ((stellar_continuum + nebular_continuum)/float(len(line)))*(3E8)/((lam/float(len(line)))**2*1E-10)

        EW = np.log10(line_luminosity/total_continuum)

        EWs.append(EW)


    if iZ==0:
        Z_label = '0.00001'
    else:
        Z_label = str(Z)

    ax.plot(grid['log10age'], np.array(EWs), label = r'$\rm Z='+Z_label+'$', c=color, ls=ls, lw=1)


ax.set_xlim([6.,9.])
ax.legend(fontsize=6)


ax.set_xlabel(r'$\rm\log_{10}(age/yr)$')
ax.set_ylabel(r'$\rm\log_{10}(EW_{'+labels[label]+'}/\AA)$')


figname = 'figures/EW_SSP_'+str(label)+'.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot
