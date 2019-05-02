

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

log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10age': log10age, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}



# ---- define plot

parameter_range = {}

parameter_range['fesc'] = [0.0, 0.1, 0.25, 0.5]
parameter_range['log10age'] = [7.,8.,9.,10.]
parameter_range['log10n_H'] = [2.,3.,4.]
parameter_range['log10U'] = [-1.,-2.,-3.,-4.]
parameter_range['CO'] = [-0.75,-0.5,-0.25,0.0,0.25]
parameter_range['d2m'] = [0.1,0.2,0.3,0.4,0.5]

for parameter in ['fesc','log10age','log10n_H', 'log10U','CO','d2m']:


    plt.style.use('simple')

    fig = plt.figure( figsize=(3,3) )

    left  = 0.2  
    bottom = 0.2   
    height = 0.7
    width = 0.7

    ax = fig.add_axes((left, bottom, width, height))


    # --- add vertical lines to denote solar metallicity and 1/5 and 1/50

    Z_sol = 0.02

    for s in [50,5,1]:

        ax.axvline(np.log10(Z_sol/float(s)), c= 'k', alpha = 0.1, lw = '1')

        if s == 1:
            ax.text(np.log10(Z_sol/float(s))-0.075, 1.0, r'$\rm Z_{\odot}$', rotation = 90., fontsize = 5., alpha =0.5, horizontalalignment='center', verticalalignment='center') 
        else:
            ax.text(np.log10(Z_sol/float(s))-0.075, 1.0, r'$\rm 1/'+str(s)+'Z_{\odot}$', rotation = 90., fontsize = 5., alpha =0.5, horizontalalignment='center', verticalalignment='center') 
            print r'$\rm 1/'+str(s)+'Z_{\odot}$'


    # --- make graph


    for i, parameter_value in enumerate(parameter_range[parameter]):

        c = cm.plasma(float(i)/len(parameter_range[parameter]))

        p = copy.copy(default_p)

        p[parameter] = parameter_value

        EW = []
        for iZ, log10Z in enumerate(log10Zs):

            p['log10Z'] = log10Z

            grid.generate_SED(p, dust_model = False)
            grid.generate_line_luminosity(p, line, dust_model = False)
            EW.append(grid.EW())
   
        # ax.semilogy(log10Zs, EW, c = c, label = r'$\rm \log_{10}(t_{SF}/yr)='+str(log10age)+'$', alpha = 0.8, lw =1.0)
        
        ax.semilogy(log10Zs, EW, c = c, label = r'$\rm'+readGrid.parameter_labels[parameter]+'='+str(parameter_value)+'$', alpha = 0.8, lw =1.0)




    # --- add model labels

    model_label = r'$\rm BPASSv2.1\quad Modified\ Salpeter\ IMF\ (m_{\star}=0.1-300\,{\rm M_{\odot}})\quad Cloudy\ (C17.00)$'

    plt.figtext(left+width/2., bottom+height+0.04, model_label, size = 5, alpha = 0.4, ha = 'center', va = 'bottom')

    # --- default parameters

    model_label = ' \quad '.join([ 'default:', r'f_{esc}='+str(fesc), r'\log_{10}(t_{SF}/yr)='+str(int(np.round(default_p['log10age'],0))), r'\log_{10}(n/cm^{-3})='+str(int(np.round(default_p['log10n_H'],0))), r'\log_{10}(U_S)='+str(int(np.round(default_p['log10U'],1))), r'\log_{10}[(C/O)/(C/O)_{\odot}]='+str(default_p['CO']), r'\xi_{d}='+str(default_p['d2m'])])
   
    plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 4, alpha = 0.4, ha = 'center', va = 'bottom')



    # --- add legend

    ax.legend(prop={'size':5}, labelspacing = 0.0)


    ax.set_xlim([-4., log10Zs[-1]])
    ax.set_ylim([0.003, 30])

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: '{:g}'.format(y))) # suppress scientific notation

    ax.set_ylabel(r'$\rm equivalent\ width/\AA$')
    ax.set_xlabel(r'$\rm\log_{10}(Z)$')
  
    fig_name = 'figures/CIII_EW_Z_'+parameter+'.pdf'
  
    print fig_name 
  
    fig.savefig(fig_name)



