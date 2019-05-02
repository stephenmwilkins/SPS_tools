


import numpy as np



lam, cID, intrinsic, emergent = np.loadtxt('lines_z020_bin_6.2_v1.lines', dtype = str, delimiter='\t', usecols = (0,1,2,3)).T


# --- Halpha luminosity

line_lam = 6563.

s = [cID=='H  1  6563A']

Ha = float(emergent[s][0])

print 'log10(Halpha)', Ha



# --- Elizabeth's method

e = np.array([cid.split(' ')[0] for cid in cID])
s = [e == 'nInu']

continuum_lam = np.array(map(float, lam[s][::-1]))
continuum_nInu = np.array(map(float, emergent[s][::-1]))

# for l, I in zip(continuum_lam,continuum_nInu): print l, I

line_continuum = np.interp(line_lam, continuum_lam, continuum_nInu)

print '------------'

print 'Elizabeths Method:', line_lam*(10**Ha)/10**line_continuum



print '------------'


lam, incident, transmitted, nebular, total, linecont = np.loadtxt('out_z020_bin_6.2_v1.cont', delimiter='\t', usecols = (0,1,2,3,4,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total

nebular_continuum = nebular - linecont

total_continuum = nebular_continuum + transmitted

for t in ['incident', 'transmitted', 'nebular_continuum', 'total_continuum']:

    globals()[t] *= 4.*np.pi * (100.*10*3.0857E16)**2 # L erg s^-1 # no longer needed in C17

    line_continuum = np.interp(line_lam, lam[::-1], globals()[t][::-1]) 

    print t, line_lam*(10**Ha)/line_continuum




import matplotlib.pyplot as plt
from matplotlib import cm

plt.style.use('simple')

fig = plt.figure( figsize=(5,4) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))

s = [(lam>1000.)&(lam<10000.)]

for t in ['incident', 'transmitted', 'nebular_continuum', 'total_continuum']:

   ax.plot(np.log10(lam[s]), np.log10(globals()[t][s]), label = t)

for l, I in zip(continuum_lam,continuum_nInu): 

    if l<10000. and l>1000.:

        ax.scatter(np.log10(l), I, c = 'k', zorder = 4, s = 4)



ax.legend(fontsize=6, labelspacing = 0.0)

ax.set_ylabel(r'$\rm\log_{10}(\nu L_{\nu})$')
ax.set_xlabel(r'$\rm\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/SED.pdf')

