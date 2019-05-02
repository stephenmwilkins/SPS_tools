

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')




SPS = 'BPASSv2.2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


# lines = [['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563']]


line_ratios = [[['OII3726','OII3729'], ['OIII5007']], [['OIII5007'],['HI4861']], [['NII6583'],['HI6563']], [['SII6731','SII6716'],['HI6563']], [['NII6583'],['OII3726','OII3729']]  ]



N_line_ratios = len(line_ratios)



fig, axes = plt.subplots(N_line_ratios, 1, figsize=(3.5, N_line_ratios*1.5), sharey = True, sharex = True)

left  = 0.1
right = 0.95
bottom = 0.2/N_line_ratios
top = 0.98
wspace = 0.02   # the amount of width reserved for blank space between subplots
hspace = 0.02

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

axes = axes.flatten()


for log10tmax in [6.,7.,8.]:

    grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS+'/'+IMF+'/lines.p','rb'), encoding='latin1')

#     for l in grid['lines']: print(l)


    for ax, line_ratio in zip(axes, line_ratios):

        print(line_ratio)

        R = []

        for iZ, Z in enumerate(grid['Z']):

            r0 = 0.0 # R = r0/r1
            r1 = 0.0 # R = r0/r1

            totSF = 0.0

            for ia, log10age in enumerate(np.arange(6.,log10tmax+0.1,0.1)):

                # --- determine SF activity in bin

                SF = 10**log10age - 10**(log10age-0.1) 
                if log10age==6.0: SF = 10**6.0     
                totSF += SF
                

                for l in line_ratio[0]:
                    l = l.encode('UTF-8')
                    r0 += SF*10**grid[l]['luminosity'][ia, iZ]

                for l in line_ratio[1]:
                    l = l.encode('UTF-8')
                    r1 += SF*10**grid[l]['luminosity'][ia, iZ]

            
            R.append(np.log10(r0/r1))

        ax.plot(grid['log10Z'], R ) #, label = r'$\rm '+labels[label]+r'$', c=c, ls=ls, lw=1)


        
        ax.set_xlim([-5.0, -1.4])
        ax.set_ylim([-4.5, 2.5])


        label = 'R='+'-'.join(line_ratio[0])+'/'+'-'.join(line_ratio[1])

        ax.text( -4.8, 1.5, label, size=7.)

        # ax.set_ylabel(r'$\rm\log_{10}('+label+')$')



axes[-1].set_xlabel(r'$\rm\log_{10}(Z)$')


figname = 'figures/line_ratios.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



