import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle

from scipy import stats

plt.style.use('simple')

fancy = lambda x: r'$\rm '+str(r'\ '.join(x.split(' ')))+' $'



lines = [[['NII6583'],['HI6563']], [['OIII5007'],['HI4861']]]

SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'
model_label = r'$\rm BPASSv2.2.1\quad Modified\,Salpeter\, (m_{\star}=0.1-300\,{\rm M_{\odot}})$'







fig = plt.figure(figsize = (3.5, 2.5))

left  = 0.15
bottom = 0.15
width = 0.55
height = 0.8

ax = fig.add_axes((left, bottom, width, height))


grid = pickle.load(open('../../BuildGrid/Z_U_nH/grids/'+SPS+'/'+IMF+'/lines_SSP1.p','rb'), encoding='latin1')

inH_def = (np.abs(grid['log10n_H'] - 2.5)).argmin()
iU_def = (np.abs(grid['log10U_S_0'] - -2.0)).argmin()


print(grid['log10U_S_0'], iU_def)
print(grid['log10n_H'], inH_def)


inHs = [(np.abs(grid['log10n_H'] - log10n_H)).argmin() for log10n_H in [1., 2.,2.5, 3., 4.]]      
iUs = [(np.abs(grid['log10U_S_0'] - log10U_S_0)).argmin() for log10U_S_0 in [-1.,-2.,-3.]]      


# R1 = []
# R2 = []
# 
# for iZ, Z in enumerate(grid['Z']):
# 
#     r0 = 0.0 # R = r0/r1
#     r1 = 0.0 # R = r0/r1
# 
#     for l in lines[0][0]:
#         l = l.encode('UTF-8')
#         r0 += grid[l]['luminosity'][iZ, iU_def, inH_def]
# 
#     for l in lines[0][1]:
#         l = l.encode('UTF-8')
#         r1 += grid[l]['luminosity'][iZ, iU_def, inH_def]
# 
#     R1.append(np.log10(r0/r1))
#     
#     r0 = 0.0 # R = r0/r1
#     r1 = 0.0 # R = r0/r1
# 
#     for l in lines[1][0]:
#         l = l.encode('UTF-8')
#         r0 += grid[l]['luminosity'][iZ, iU_def, inH_def]
# 
#     for l in lines[1][1]:
#         l = l.encode('UTF-8')
#         r1 += grid[l]['luminosity'][iZ, iU_def, inH_def]
# 
#     R2.append(np.log10(r0/r1))


for iU, ls in zip(iUs,[':','-','--']):

    R1 = []
    R2 = []
    
    for iZ, Z in enumerate(grid['Z']):

        r0 = 0.0 # R = r0/r1
        r1 = 0.0 # R = r0/r1

        for l in lines[0][0]:
            l = l.encode('UTF-8')
            r0 += grid[l]['luminosity'][iZ, iU, inH_def]

        for l in lines[0][1]:
            l = l.encode('UTF-8')
            r1 += grid[l]['luminosity'][iZ, iU, inH_def]

        R1.append(np.log10(r0/r1))
    
        r0 = 0.0 # R = r0/r1
        r1 = 0.0 # R = r0/r1

        for l in lines[1][0]:
            l = l.encode('UTF-8')
            r0 += grid[l]['luminosity'][iZ, iU, inH_def]

        for l in lines[1][1]:
            l = l.encode('UTF-8')
            r1 += grid[l]['luminosity'][iZ, iU, inH_def]

        R2.append(np.log10(r0/r1))


    ax.plot(R1, R2 , c='k', lw=1, ls=ls, zorder = 0) #, label = r'$\rm '+labels[label]+r'$', c=c, ls=ls, lw=1)
    ax.scatter(R1, R2, c='k', s=10, zorder = 1)
    
    c = cm.viridis(np.arange(len(grid['Z']))/len(grid['Z']))
    
    ax.scatter(R1, R2, c=c, s=5, zorder = 2)




# --- SDSS


from astropy.io import fits

SDSS_data = fits.open('obs_data/SDSS_line_fluxes_BPT.fits')[1]


sel = (SDSS_data.data['nii_6584_flux']>0) & (SDSS_data.data['h_alpha_flux']>0) & \
      (SDSS_data.data['oiii_5007_flux']>0) & (SDSS_data.data['h_beta_flux']>0)


NII_HALPHA =  np.log10( SDSS_data.data['nii_6584_flux'][sel] / SDSS_data.data['h_alpha_flux'][sel] )
OIII_HBETA =  np.log10( SDSS_data.data['oiii_5007_flux'][sel] / SDSS_data.data['h_beta_flux'][sel] )

gridsize = 100
extent = (-3.,1.,-2.,2.)


ax.hexbin(NII_HALPHA, OIII_HBETA, gridsize = gridsize, bins = 'log', cmap='Greys', extent = extent, linewidths=0., mincnt = 1, zorder = 0, alpha = 0.5)



def Kewley2001(NII):
    OIII = (0.61/( NII-0.47 )) + 1.19
    OIII[NII > 0.47] = -100
    return OIII

def Kauffmann2003(NII):
    OIII = (0.61/(NII-0.05))+1.3
    OIII[NII > 0.05] = -100
    return OIII

NII = np.linspace(-4,2,1000)

# ax.plot( NII ,Kewley2001(NII)      ,color='dodgerblue', lw = 1, linestyle='--' ,label='Kewley+ 2001', zorder=1, alpha = 0.5)
# ax.plot( NII ,Kauffmann2003(NII)   ,color='dodgerblue' , lw = 1, label='Kauffmann+ 2003' , zorder=1, alpha = 0.5)





ax.set_xlim([-4.,1.])
ax.set_ylim([-2.,2.])

ax.set_xlabel(r'$\rm \log_{10}(\regular{[NII]6583}/H\alpha)$')
ax.set_ylabel(r'$\rm \log_{10}(\regular{[OIII]5007}/H\beta)$')




dummies = []
labels = []
# dummies += [ax.plot([], [], ls='-', c='k', alpha=0.5, lw = lw)[0] for lw in [1,2]]        

labels += [r'$\rm\log_{10}(U_{S,ref})='+str(x)+'$' for x in [-1.,-2.,-3.]]
dummies += [ax.plot([], [], ls=ls, c='k', alpha=0.5, lw = 1)[0] for ls in [':','-','--']]        

labels += [r'$\rmZ='+str(x)+'$' for x in grid['Z']]
dummies += [ax.scatter([], [], c=c, s=10, lw=0.5, edgecolor='k') for c in cm.viridis(np.arange(len(grid['Z']))/len(grid['Z']))]        


ax.legend(dummies, labels,loc = 'center left', bbox_to_anchor = (1.0, 0.5), fontsize = 6, handlelength=2.5)










fig.savefig('figures/BPT_theory.pdf')

fig.clf()