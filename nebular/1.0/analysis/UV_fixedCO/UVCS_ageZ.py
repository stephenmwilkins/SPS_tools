

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import matplotlib.lines as mlines

import readGrid
import copy



SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'


grid = readGrid.grid(SPS, IMF, lines = False, generate_SED = True)




Zs = np.array([0.00001, 0.0001, 0.001, 0.004, 0.01, 0.04])

log10ages = np.arange(7.,10.,0.1)

# ---- define default parameters

fesc = 0.0

log10age = 8.
log10Z = np.log10(0.004)

log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10Z': log10Z, 'log10age': log10age, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}



plt.style.use('simple')

fig = plt.figure( figsize=(3,3) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))




for i, Z in enumerate(Zs):

    UVCS = []

    c = cm.plasma(float(i)/len(Zs))

    for log10age in log10ages:

        p = copy.copy(default_p)

        p['log10age'] = log10age
        p['log10Z'] = np.log10(Z)
  
        grid.generate_SED(p, dust_model = False)
        grid.generate_UVCS()
        UVCS.append(grid.UVCS_int)


    precision = int(np.ceil(abs(np.log10(abs(Z))))) if p != 0 else 1

    ax.plot(log10ages, UVCS, label = '$'+'{:.{}f}'.format(Z, precision)+'$', c = c, alpha = 0.8, lw =1.0)



for i, Z in enumerate(Zs):

    UVCS = []

    c = cm.plasma(float(i)/len(Zs))

    for log10age in log10ages:

        p = copy.copy(default_p)

        p['fesc'] = 1.0
        p['log10age'] = log10age
        p['log10Z'] = np.log10(Z)
  
        grid.generate_SED(p, dust_model = False)
        grid.generate_UVCS()
        UVCS.append(grid.UVCS_int)

    ax.plot(log10ages, UVCS, c = c, alpha = 0.8, lw =1.0, ls = '-.')









# --- add model labels

model_label = r'$\rm BPASSv2.1\quad Modified\ Salpeter\ IMF\ (m_{\star}=0.1-300\,{\rm M_{\odot}})\quad Cloudy\ (C17.00)$'

plt.figtext(left+width/2., bottom+height+0.04, model_label, size = 5, alpha = 0.4, ha = 'center', va = 'bottom')

# --- default parameters

model_label = ' \quad '.join([ 'default:', r'f_{esc}='+str(fesc), r'\log_{10}(n/cm^{-3})='+str(int(np.round(default_p['log10n_H'],0))), r'\log_{10}(U_S)='+str(int(np.round(default_p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(default_p['CO']), r'\xi_{d}='+str(default_p['d2m'])])

plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 4, alpha = 0.4, ha = 'center', va = 'bottom')





# --- add escape fraction legend

fesc0 = mlines.Line2D([], [], lw = 1, color='0.5', ls = '-', label=r'$\rm f_{esc}=0$')
fesc1 = mlines.Line2D([], [], lw = 1, color='0.5', ls = '-.', label=r'$\rm f_{esc}=1$')

fesc_legend = ax.legend(loc = 'lower right', handles=[fesc0, fesc1], handlelength = 2.5, fontsize=6, labelspacing = 0.0)
leg = plt.gca().add_artist(fesc_legend)

# --- add legend

ax.legend(prop={'size':5}, labelspacing = 0.0, loc='upper left')


ax.set_xlim([7.,9.5])
ax.set_ylim([-3.5, -1.5])

ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y))) # suppress scientific notation

ax.set_ylabel(r'$\rm UV\ continuum\ slope\ (\beta)$')
ax.set_xlabel(r'$\rm \log_{10}(t_{SF}/yr)$')










fig_name = 'figures/UVCS_ageZ.pdf'







print fig_name 

fig.savefig(fig_name)



