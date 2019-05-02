

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker

import readGrid
import copy



SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'


line = 'CIII1907,CIII1909'
line_label = 'CIII]'

grid = readGrid.grid(SPS, IMF, [line], generate_SED = True)




# Zs = np.array([0.00001, 0.0001, 0.001, 0.004, 0.01, 0.04])
# log10Zs = np.log10(Zs)

log10ages = grid.stellar_grid['log10age']

log10Zs = grid.stellar_grid['log10Z']

print log10Zs


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





# --- CF00

tau_V_hats = np.arange(0.0, 4.0, 0.01)

mus = [0.0,0.25,0.5,0.75,1.0]

for i, mu in enumerate(mus):

    EW = []
    A1500 = []

    print mu

    c = cm.plasma(float(i)/len(mus))

    for tau_V_hat in tau_V_hats:

        p = copy.copy(default_p)

        p['tau_V_hat'] = tau_V_hat
        p['mu'] = mu
  
        grid.generate_SED(p, dust_model = 'CF00')
        grid.generate_line_luminosity(p, line, dust_model = 'CF00')  
        EW.append(grid.EW())
        A1500.append(grid.A1500())

        # print tau_V_hat, grid.A1500(), grid.EW()

    ax.semilogy(A1500, EW, c = c, label = r'$\rm\mu='+str(mu)+'$', alpha = 0.8, lw =1.0)


# --- Calzetti

tau_V_nebulars = np.arange(0.0, 4.0, 0.01)

EW = []
A1500 = []

for tau_V_nebular in tau_V_nebulars:

    p = copy.copy(default_p)

    p['tau_V_nebular'] = tau_V_nebular
    

    grid.generate_SED(p, dust_model = 'Calzetti')
    grid.generate_line_luminosity(p, line, dust_model = 'Calzetti')  
    EW.append(grid.EW())
    A1500.append(grid.A1500())

ax.semilogy(A1500, EW, c = 'k', label = r'$\rm Calzetti\ et\ al.\ (2000)$', alpha = 0.4, lw =1.5)



# --- add model labels

model_label = r'$\rm BPASSv2.1\quad Modified\ Salpeter\ IMF\ (m_{\star}=0.1-300\,{\rm M_{\odot}})\quad Cloudy\ (C17.00)$'

plt.figtext(left+width/2., bottom+height+0.04, model_label, size = 5, alpha = 0.4, ha = 'center', va = 'bottom')

# --- default parameters

model_label = ' \quad '.join([ 'default:', r'f_{esc}='+str(fesc), r'\log_{10}(t_{SF}/yr)='+str(int(np.round(default_p['log10age'],0))), r'\log_{10}(n/cm^{-3})='+str(int(np.round(default_p['log10n_H'],0))), r'\log_{10}(U_S)='+str(int(np.round(default_p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(default_p['CO']), r'\xi_{d}='+str(default_p['d2m'])])

plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 4, alpha = 0.4, ha = 'center', va = 'bottom')



# --- add legend

ax.legend(prop={'size':5}, labelspacing = 0.0)


ax.set_xlim([0.,2.5])
ax.set_ylim([0.1, 20])

ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y))) # suppress scientific notation

ax.set_ylabel(r'$\rm equivalent\ width/\AA$')
ax.set_xlabel(r'$\rm A_{1500}$')

fig_name = 'figures/CIII_EW_dust.pdf'

print fig_name 

fig.savefig(fig_name)



