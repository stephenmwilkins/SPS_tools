

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




# ---- define default parameters

fesc = 0.0

log10age = 8.
log10Z = np.log10(0.004)

log10n_H = 2.0
log10U = -2.0
d2m = 0.3
CO = 0.0

default_p = {'fesc': fesc, 'log10Z': log10Z, 'log10age': log10age, 'log10n_H': log10n_H, 'log10U': log10U, 'd2m': d2m, 'CO':CO}


p = copy.copy(default_p)

grid.generate_nebular_continuum(p)
grid.generate_stellar(p)
grid.generate_line_luminosity(p, line)  

print 'EW:', grid.EW()


p = copy.copy(default_p)

p['tau_V_hat'] = 0.1
p['mu'] = 0.5

grid.generate_nebular_continuum(p)
grid.generate_stellar(p)
grid.generate_line_luminosity(p, line)  


print grid.stellar/grid.intrinsic_stellar
# print grid.nebular_continuum/grid.intrinsic_nebular_continuum

print 'fesc(1500):', grid.fesc(1500.)
print 'A1500:', grid.A1500()
print 'EW:', grid.EW()

  
