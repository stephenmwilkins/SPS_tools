

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



ax.fill_between([-3.0,-2.0], [-7,-7], [5.,5.], color='k',alpha=0.05)



SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'



lines = [['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563']]

# lines = [['HI6563']]


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


for log10tmax in [7.]:

    grid = pickle.load(open('/Users/stephenwilkins/research/data/SPS/nebular/2.0/Z_refQ_wdust/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')

    for il, line, ls in zip(range(len(lines)), lines, ['-','-.','--',':']*3):

        c = cm.magma(float(il)/float(len(lines)))

        EWs = []

        label = '-'.join(line)

        for iZ, Z in enumerate(grid['Z']):

            line_luminosity = 0.0
            continuum = 0.0
            lam = 0.0

            totSF = 0.0

            for ia, log10age in enumerate(np.arange(6.,log10tmax+0.1,0.1)):

                # --- determine SF activity in bin

                SF = 10**log10age - 10**(log10age-0.1)
                if log10age==6.0: SF = 10**6.0
                totSF += SF

                print(log10age, SF, totSF)

                for l in line:

                    # l = l.encode('UTF-8')

                    line_luminosity += SF*10**grid[l]['luminosity'][ia, iZ]
                    continuum += SF*grid[l]['total_continuum'][ia, iZ]

                    if ia == 0: lam += grid[l]['lam']


            total_continuum = ((continuum)/float(len(line)))*(3E8)/((lam/float(len(line)))**2*1E-10)

            EW = np.log10(line_luminosity/total_continuum)

            print(Z, EW)

            EWs.append(EW)

        ax.plot(grid['log10Z'], np.array(EWs), label = r'$\rm '+labels[label]+r'$', c=c, ls=ls, lw=1)



# ax.legend(loc = 'center left', fontsize=7, bbox_to_anchor = (1.0, 0.5))

ax.legend(loc = 'lower left', fontsize=6, handlelength=2.5, labelspacing=0.2)


ax.set_xlim([-5.0, -1.4])
ax.set_ylim([-3.5, 3.6])


ax.set_xlabel(r'$\rm\log_{10}(Z)$')
ax.set_ylabel(r'$\rm\log_{10}(EW/\AA)$')


figname = 'figures/lines_EWs.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot
