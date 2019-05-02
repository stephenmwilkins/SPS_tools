

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker

import readGrid
import copy


SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'



lines = ['CIV1548,CIV1551','CIII1907,CIII1909']
line_labels = [r'C\regular_{IV}',r'[C\regular_{III}],C\regular_{III}]']

grid = readGrid.grid(SPS, IMF, lines, generate_SED = False)


log10Zs = grid.nebular_line_grid['log10Z']
    
# ---- define default parameters

fesc = 0.0
log10n_H = 1.0
log10U = -1.5
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}



# ---- define plot


plt.style.use('simple')


npanels = 4

fig, axes = plt.subplots(1, npanels, figsize=(6.75, 2.25), sharey = True, sharex = True)

left  = 0.1
right = 0.95
bottom = 0.2
top = 0.8
wspace = 0.02   # the amount of width reserved for blank space between subplots

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace)





for ax, parameter in zip(axes.flatten(), ['log10U', 'log10n_H', 'd2m', 'CO']):

    for iv, v in enumerate(grid.nebular_line_grid[parameter]):

        p = copy.copy(default_p)
    
        p[parameter] = v

        c = cm.plasma(0.1+float(iv)/len(grid.nebular_line_grid[parameter]))

        R = []
        for iZ, log10Z in enumerate(log10Zs):

            p['log10Z'] = log10Z

            grid.generate_line_luminosity(p, lines[0], dust_model = False)

            L0 = grid.line_luminosity

            grid.generate_line_luminosity(p, lines[1], dust_model = False)

            L1 = grid.line_luminosity

            R.append(np.log10(L0/L1))


        ax.plot(log10Zs, R, c = c, label = '$'+str(v)+'$', alpha = 0.8, lw =1.0)


    leg = ax.legend(loc='lower left', scatterpoints = 1, prop={'size':4}, labelspacing = 0.0)

    ax.text(0.5, 1.01, r'$\rm'+readGrid.parameter_labels[parameter]+'$', fontsize=6, ha = 'center', va='bottom', transform=ax.transAxes)

    ax.set_xlim([-4., log10Zs[-1]])
    ax.set_ylim([-2,2.])


  
axes[0].set_ylabel(r'$\rm \log_{10}('+'/'.join(line_labels)+')$', fontsize=7)
fig.text(0.525, 0.05, r'$\rm\log_{10}(Z)$', ha='center', va='center', fontsize=7)

fig.text(0.525, 0.9, model_label, ha='center', va='bottom', fontsize=6, alpha=0.4)  
  
fig_name = 'figures/line_ratios.pdf'
  
print fig_name 
  
fig.savefig(fig_name)



