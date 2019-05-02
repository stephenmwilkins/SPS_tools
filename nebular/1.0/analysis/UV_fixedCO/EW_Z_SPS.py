

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import matplotlib.lines as mlines

import readGrid
import copy


run = True


SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'


line = 'CIII1907,CIII1909'
line_label = 'CIII]'


    
# ---- define default parameters

fesc = 0.0

log10age = 8.

log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10age': log10age, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}



# ---- define plot


plt.style.use('simple')

fig = plt.figure( figsize=(3,3) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))


Z_sol = 0.02

for s in [50,5,1]:

    ax.axvline(np.log10(Z_sol/float(s)), c= 'k', alpha = 0.1, lw = '1')

    
    if s == 1:
        ax.text(np.log10(Z_sol/float(s))-0.075, 1.0, r'$\rm Z_{\odot}$', rotation = 90., fontsize = 5., alpha =0.5, horizontalalignment='center', verticalalignment='center') 
    else:
        ax.text(np.log10(Z_sol/float(s))-0.075, 1.0, r'$\rm 1/'+str(s)+'Z_{\odot}$', rotation = 90., fontsize = 5., alpha =0.5, horizontalalignment='center', verticalalignment='center') 
        print r'$\rm 1/'+str(s)+'Z_{\odot}$'



SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = True)
log10Zs = grid.stellar_grid['log10Z']

EW = []
for iZ, log10Z in enumerate(log10Zs):

    p = copy.copy(default_p)

    p['log10Z'] = log10Z
    
    grid.generate_SED(p, dust_model = False)
    grid.generate_line_luminosity(p, line, dust_model = False)
    EW.append(grid.EW())

ax.semilogy(log10Zs, EW, c = 'k', alpha = 0.6, lw =1.0, ls = '-')




SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_100'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = True)
log10Zs = grid.stellar_grid['log10Z']

EW = []
for iZ, log10Z in enumerate(log10Zs):

    p = copy.copy(default_p)

    p['log10Z'] = log10Z
    
    grid.generate_SED(p, dust_model = False)
    grid.generate_line_luminosity(p, line, dust_model = False)
    EW.append(grid.EW())

ax.semilogy(log10Zs, EW, c = 'k', alpha = 0.6, lw =1.0, ls = '--')




SPS = 'P2'
IMF = 'ModSalpeter_100'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = True)
log10Zs = grid.stellar_grid['log10Z']

EW = []
for iZ, log10Z in enumerate(log10Zs):

    p = copy.copy(default_p)

    p['log10Z'] = log10Z
    
    grid.generate_SED(p, dust_model = False)
    grid.generate_line_luminosity(p, line, dust_model = False)
    EW.append(grid.EW())

ax.semilogy(log10Zs, EW, c = 'k', alpha = 0.6, lw =1.0, ls = '-.')



# --- add SPS mode legend

mod1 = mlines.Line2D([], [], lw = 1, color='0.5', ls = '-', label=r'$\rm BPASSv2.1\ Modified\ Salpeter\ (m=0.1-300\, M_{\odot})$')
mod2 = mlines.Line2D([], [], lw = 1, color='0.5', ls = '--', label=r'$\rm BPASSv2.1\ Modified\ Salpeter\ (m=0.1-100\, M_{\odot})$')
mod3 = mlines.Line2D([], [], lw = 1, color='0.5', ls = '-.', label=r'$\rm PEGASE.2\ Modified\ Salpeter\ (m=0.1-100\, M_{\odot})$')

fesc_legend = ax.legend(loc = 'lower left', handles=[mod1, mod2, mod3], handlelength = 2.5, fontsize=5, labelspacing = 0.0)
leg = plt.gca().add_artist(fesc_legend)







model_label = ' \quad '.join([ 'default:', r'f_{esc}='+str(fesc), r'\log_{10}(t_{SF}/yr)='+str(int(np.round(default_p['log10age'],0))), r'\log_{10}(n/cm^{-3})='+str(int(np.round(default_p['log10n_H'],0))), r'\log_{10}(U_S)='+str(int(np.round(default_p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(default_p['CO']), r'\xi_{d}='+str(default_p['d2m'])])
   
plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 4, alpha = 0.4, ha = 'center', va = 'bottom')

# ax.legend(prop={'size':6}, labelspacing = 0.0, title = r'$\rm\log_{10}(Z)$')
# ax.legend(prop={'size':6}, labelspacing = 0.0)



ax.set_xlim([-4., log10Zs[-1]])
ax.set_ylim([0.003, 30])

ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y))) # suppress scientific notation

ax.set_ylabel(r'$\rm equivalent\ width/\AA$')
ax.set_xlabel(r'$\rm\log_{10}(Z)$')
  
fig_name = 'figures/CIII_EW_Z_SPS.pdf'
  
print fig_name 
  
fig.savefig(fig_name)



