

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import readGrid



SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

lines = ['HI6563']

grid = readGrid.grid(SPS, IMF, lines, generate_SED = True)


fesc = 0.0
log10age = 8.05 # --- must be *.*5
# log10age = 7.05 # --- must be *.*5
Z = 0.0001
log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

p = {'fesc': fesc, 'log10age': log10age, 'log10Z':np.log10(Z), 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}


print grid.line_luminosity(p, 'HI6563')

'KE12=41.27'

grid.generate_nebular_continuum(p)

grid.generate_stellar(p)

print np.log10(grid.EW(p, 'HI6563'))




plt.style.use('simple')

fig = plt.figure( figsize=(6,4) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))


s = [(grid.lam<10000)&(grid.lam>1000)]

ax.plot(np.log10(grid.lam[s]), np.log10(grid.intrinsic_stellar[s]), label = 'intrinsic stellar', c = 'k', alpha = 0.8, lw =0.5)

ax.plot(np.log10(grid.lam[s]), np.log10(grid.intrinsic_nebular_continuum[s]), label = 'intrinsic nebular', c = 'k', alpha = 0.8, lw =0.5)


ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/test_SED.pdf')