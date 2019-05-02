

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import readGrid
import copy





SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'


line = 'HI6563'
line = 'OII3729,OII3726'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = False)

log10Zs = grid.nebular_line_grid['log10Z']

log10Zs = [-5.,-4.,-3.,-2.69897,-2.22184875,-2.09691001,-2.,-1.69897,-1.52287875,-1.39794001]

print log10Zs


# ---- define default parameters

fesc = 0.0
Z = 0.01
log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10Z':np.log10(Z), 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}

default_luminosity = grid.line_luminosity(default_p, line)


# ---- define plot


plt.style.use('simple')



npanels = 2

fig, axes = plt.subplots(1, npanels, figsize=(3, 1.5), sharey = True, sharex = True)

left  = 0.1
right = 0.95
bottom = 0.2
top = 0.8
wspace = 0.02   # the amount of width reserved for blank space between subplots

fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace)





# ------ 


for ax, parameter in zip(axes.flatten(), ['log10U', 'log10n_H']):

    p = copy.copy(default_p)

    for iv, v in enumerate(grid.nebular_line_grid[parameter]):

        c = cm.plasma(0.1+float(iv)/len(grid.nebular_line_grid[parameter]))

        p[parameter] = v

        rluminosity = []

        for log10Z in log10Zs:

            p['log10Z'] = log10Z

            luminosity = grid.line_luminosity(p, line)

            rluminosity.append(np.log10(luminosity/default_luminosity))
            
        ax.plot(log10Zs, rluminosity, label = '$'+str(v)+'$', c = c, alpha = 0.8, lw =0.5)
  
  
    leg = ax.legend(loc='lower right', scatterpoints = 1, prop={'size':4}, labelspacing = 0.0)

    ax.text(0.5, 1.01, r'$\rm'+readGrid.parameter_labels[parameter]+'$', fontsize=6, ha = 'center', va='bottom', transform=ax.transAxes)


axes[0].set_ylabel(r'$\rm\log_{10}(L/L_{default})$', fontsize=7)
fig.text(0.525, 0.05, r'$\rm\log_{10}(Z)$', ha='center', va='center', fontsize=7)


fig.text(0.525, 0.9, model_label, ha='center', va='bottom', fontsize=6, alpha=0.4)

#     ax.set_ylabel(r'$L_{\nu}$')
#     ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/line_luminosities_'+str(line)+'.pdf')