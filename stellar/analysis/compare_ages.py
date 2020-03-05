

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
import cPickle as pickle


# SPS = 'P2'
# IMF = 'ModSalpeter_100'

SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'


plt.style.use('simple')

fig = plt.figure( figsize=(6,4) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))



stellar_grid = pickle.load(open('../../stellar/BuildGrid/grids/'+SPS+'/'+IMF+'/stellar.p','r'))

print stellar_grid.keys()

print stellar_grid['L_nu'].shape


Z = 0.008

print stellar_grid['log10age']
print stellar_grid['Z']

iZ = np.where(stellar_grid['Z']==Z)[0][0]

print iZ

lam = stellar_grid['lam']


ages = np.arange(6.0, 10.0, 0.5)

for i,a in enumerate(ages):

    print a, np.where(stellar_grid['log10age']==a)

    c = cm.plasma(float(i)/float(len(ages)+1))
    
    # ia = np.where(stellar_grid['log10age']==a)[0][0]
    
    ia = (np.fabs(stellar_grid['log10age']-a)).argmin()
    
    print ia
    
    Lnu = stellar_grid['L_nu'][ia, iZ]

    s = [(lam<10000)&(lam>250)]

    ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = a, c = c, lw = 1)
    










ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/ages.pdf')