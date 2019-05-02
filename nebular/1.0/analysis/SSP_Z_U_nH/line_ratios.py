

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')




SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


# lines = [['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563']]


line_ratios = [[['OII3726'], ['OIII5007']], [['OIII5007'],['HI4861']], [['NII6583'],['HI6563']], [['SII6731','SII6716'],['HI6563']], [['NII6583'],['OII3726']]  ]



labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'\regular{[CII]1907,1909}'
labels['OII3726'] = r'\regular{[OII]3727}'
labels['NeIII3869'] = r'\regular{[NeIII]3869}'
labels['NeIII3967'] = r'\regular{[NeIII]3967}'
labels['HI4340'] = r'\regular{H\gamma}'
labels['HI4861'] = r'\regular{H\beta}'
labels['OIII5007'] = r'\regular{[OIII]5007}'
labels['HI6563'] = r'\regular{H\alpha}'
labels['NII6583'] = r'\regular{[NII]6584}'
labels['SII6731-SII6716'] = r'\regular{[SII]6717,6731}'

N_line_ratios = len(line_ratios)



fig, axes = plt.subplots(1, N_line_ratios, figsize=(N_line_ratios*1.5, 2), sharey = True, sharex = True)

left  = 0.3/N_line_ratios
right = 0.85
bottom = 0.15
top = 0.9
wspace = 0.02   # the amount of width reserved for blank space between subplots
hspace = 0.02

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

axes = axes.flatten()


grid = pickle.load(open('../../BuildGrid/Z_U_nH/grids/'+SPS+'/'+IMF+'/lines_SSP1.p','rb'), encoding='latin1')

inH_def = (np.abs(grid['log10n_H'] - 2.5)).argmin()
iU_def = (np.abs(grid['log10U_S_0'] - -2.0)).argmin()


print(grid['log10U_S_0'], iU_def)
print(grid['log10n_H'], inH_def)


inHs = [(np.abs(grid['log10n_H'] - log10n_H)).argmin() for log10n_H in [1., 2.,2.5, 3., 4.]]      
iUs = [(np.abs(grid['log10U_S_0'] - log10U_S_0)).argmin() for log10U_S_0 in [-1.,-2.,-3.]]      



for i, ax, line_ratio in zip(range(N_line_ratios), axes, line_ratios):


    

    print(line_ratio)

    c = cm.ocean(i/N_line_ratios)


    # ----- default model

    R = []

    for iZ, Z in enumerate(grid['Z']):

        r0 = 0.0 # R = r0/r1
        r1 = 0.0 # R = r0/r1

        for l in line_ratio[0]:
            l = l.encode('UTF-8')
            r0 += grid[l]['luminosity'][iZ, iU_def, inH_def]

        for l in line_ratio[1]:
            l = l.encode('UTF-8')
            r1 += grid[l]['luminosity'][iZ, iU_def, inH_def]
 
        R.append(np.log10(r0/r1))

    ax.plot(grid['log10Z'], R , c=c, lw=1 ) #, label = r'$\rm '+labels[label]+r'$', c=c, ls=ls, lw=1)

    
    
    for iU, ls in zip(iUs,[':','-','--']):
    
        R = []
    
        for iZ, Z in enumerate(grid['Z']):

            r0 = 0.0 # R = r0/r1
            r1 = 0.0 # R = r0/r1

            for l in line_ratio[0]:
                l = l.encode('UTF-8')
                r0 += grid[l]['luminosity'][iZ, iU, inH_def]

            for l in line_ratio[1]:
                l = l.encode('UTF-8')
                r1 += grid[l]['luminosity'][iZ, iU, inH_def]

            R.append(np.log10(r0/r1))

        ax.plot(grid['log10Z'], R, ls = ls, c=c, lw=1 ) 


    
    ax.set_xlim([-5.0, -1.4])
    ax.set_ylim([-4.5, 2.5])

    ax.fill_between([-3.0,-2.0], [-10,-10], [10,10], color='k',alpha=0.025)


    label = r'$'+labels['-'.join(line_ratio[0])]+'/'+labels['-'.join(line_ratio[1])]+'$'

    ax.text( -3.2, 2.6, label, size=7., va = 'bottom', ha='center', color=c)

    # ax.set_ylabel(r'$\rm\log_{10}('+label+')$')




dummies = []
labels = []
# dummies += [ax.plot([], [], ls='-', c='k', alpha=0.5, lw = lw)[0] for lw in [1,2]]        

labels += [r'$\rm\log_{10}(U_{S,0})='+str(x)+'$' for x in [-1.,-2.,-3.]]
dummies += [ax.plot([], [], ls=ls, c='k', alpha=0.5, lw = 1)[0] for ls in [':','-','--']]        

axes[-1].legend(dummies, labels,loc = 'center left', bbox_to_anchor = (1.0, 0.5), fontsize = 6, handlelength=2.5)



fig.text(0.5*(right - left) + left, 0.01, r'$\rm\log_{10}(Z)$', ha = 'center', va = 'bottom')
fig.text(left/2., 0.5*(top-bottom)+bottom, r'$\rm\log_{10}(Z)$', ha = 'center', va = 'center', rotation=90.)





figname = 'figures/line_ratios.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



