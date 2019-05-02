

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as mticker

import readGrid
import copy


SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'


# common lines

# 'CIII1907,CIII1909'
# 'HeII1640'
# 'NV1239,NV1243'
# 'OIII1666,OIII1661'
# 'CIV1548,CIV1551'
# 'SiIII1883,SiIII1892'

# NIII1750/HeII # nope
# SiIII1888/HeII

line_labels = {}

line_labels['CIII1907,CIII1909'] = r'C\regular_{III}'
line_labels['HeII1640'] = r'H\regular_{II}'
line_labels['NV1239,NV1243'] = r'N\regular_{V}'
line_labels['OIII1666,OIII1661'] = r'O\regular_{III}'
line_labels['CIV1548,CIV1551'] = r'C\regular_{IV}'
line_labels['SiIII1883,SiIII1892'] = r'Si\regular_{III}'



ratio_lines = [['OIII1666,OIII1661','HeII1640'],['CIII1907,CIII1909','HeII1640']]
# ratio_lines = [['SiIII1883,SiIII1892','HeII1640'],['CIII1907,CIII1909','HeII1640']]

lines = [item for sublist in ratio_lines for item in sublist] # flatten list of lines






grid = readGrid.grid(SPS, IMF, lines, generate_SED = False)


    
# ---- define default parameters

fesc = 0.0
log10n_H = 3.0
d2m = 0.3
CO = 0.0




# ---- define plot

plt.style.use('simple')

fig = plt.figure( figsize=(3,3) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))



log10Zs = [-4., -3.5, -3., -2.5, -2.]

for i,log10Z in enumerate(log10Zs):

    c = cm.plasma(0.1+float(i)/len(log10Zs))

    default_p = {'fesc': fesc, 'log10Z': log10Z, 'log10n_H': log10n_H,  'd2m': d2m, 'CO':CO}  

    X = []
    Y = []
    
    for log10U in grid.nebular_line_grid['log10U']:


        p = copy.copy(default_p)

        p['log10U'] = log10U

        grid.generate_line_luminosity(p, ratio_lines[0][0], dust_model = False)
        X0 = grid.line_luminosity
    
        grid.generate_line_luminosity(p, ratio_lines[0][1], dust_model = False)
        X1 = grid.line_luminosity
    
        grid.generate_line_luminosity(p, ratio_lines[1][0], dust_model = False)
        Y0 = grid.line_luminosity
    
        grid.generate_line_luminosity(p, ratio_lines[1][1], dust_model = False)
        Y1 = grid.line_luminosity

        X.append(np.log10(X0/X1))
        Y.append(np.log10(Y0/Y1))

    ax.plot(X, Y, c = c, alpha = 0.5, lw =1.0, ls = '-')






ax.set_xlabel(r'$\rm \log_{10}('+'/'.join([line_labels[ratio_lines[0][0]],line_labels[ratio_lines[0][1]]])+')$', fontsize=7)
ax.set_ylabel(r'$\rm \log_{10}('+'/'.join([line_labels[ratio_lines[1][0]],line_labels[ratio_lines[1][1]]])+')$', fontsize=7)


ax.set_xlim([-3.,1.5])
ax.set_ylim([-2.5,2.5])
  
fig_name = 'figures/diagnostic.pdf'
  
print fig_name 
  
fig.savefig(fig_name)



