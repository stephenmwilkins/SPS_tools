import numpy as np
import matplotlib.pyplot as plt
import os
import utilities as uu
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colorbar
from astropy.io import fits

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = 'computer modern roman'
plt.rcParams["font.size"]   = 12
plt.rcParams["font.weight"] = 700
plt.rcParams["font.weight"] = 700




def Kewley2001(NII):
    OIII = (0.61/( NII-0.47 )) + 1.19
    OIII[NII > 0.47] = -100
    return OIII

def Kauffmann2003(NII):
    OIII = (0.61/(NII-0.05))+1.3
    OIII[NII > 0.05] = -100
    return OIII

NII = np.linspace(-4,2,1000)



SDSS_data = fits.open('SDSS_line_fluxes_BPT.fits')[1]


sel = (SDSS_data.data['nii_6584_flux']>0) & (SDSS_data.data['h_alpha_flux']>0) & \
      (SDSS_data.data['oiii_5007_flux']>0) & (SDSS_data.data['h_beta_flux']>0)


NII_HALPHA =  np.log10( SDSS_data.data['nii_6584_flux'][sel] / SDSS_data.data['h_alpha_flux'][sel] )
OIII_HBETA =  np.log10( SDSS_data.data['oiii_5007_flux'][sel] / SDSS_data.data['h_beta_flux'][sel] )


fig, ax1 = plt.subplots(figsize=(7.5,6))

plt.xlim(-3.01,1.01)
plt.ylim(-2.01,2.01)


plt.scatter(  NII_HALPHA, OIII_HBETA , s=1, color='k',alpha=0.01 )
plt.plot( NII ,Kewley2001(NII)      ,color='r',linestyle='--' ,label='Kewley+ 2001' , linewidth=2.0 ,zorder=15)
plt.plot( NII ,Kauffmann2003(NII)   ,color='r' ,label='Kauffmann+ 2003' ,zorder=16 )

plt.legend(loc='best',fancybox=True,framealpha=0.5)

plt.savefig('./BPT.png')
plt.show()
