

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import pickle

plt.style.use('simple')



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15  
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))



ax.fill_between([-3.0,-2.0], [-7,-7], [5.,5.], color='k',alpha=0.05)



SPS_models = ['BPASSv2.2.1.binary/ModSalpeter_300', 'BPASSv2.2.1.binary/ModSalpeter_100','BPASSv2.2.1.binary/1p0_300','BPASSv2.2.1.binary/1p7_300']
line_styles = ['-','-.',':','--']

SPS_labels = {}


SPS_labels['BPASSv2.2.1.binary/ModSalpeter_300'] = r'$\rm BPASSv2.2.1\ Modified\ Salpeter\ IMF\ m_{up}=300\ M_{\odot}$'
SPS_labels['BPASSv2.2.1.binary/ModSalpeter_100'] = r'$\rm BPASSv2.2.1\ Modified\ Salpeter\ IMF\ m_{up}=100\ M_{\odot}$'
SPS_labels['BPASSv2.2.1.binary/1p0_300'] = r'$\rm BPASSv2.2.1\ \alpha=2.0\ m_{up}=300\ M_{\odot}$'
SPS_labels['BPASSv2.2.1.binary/1p7_300'] = r'$\rm BPASSv2.2.1\ \alpha=2.7\ m_{up}=300\ M_{\odot}$'
SPS_labels['BPASSv2.2.1.binary/Chabrier_300'] = r'$\rm BPASSv2.2.1\ Chabrier\ IMF\ m_{up}=300\ M_{\odot}$'



lines = [['HI1216'],['CIII1909','CIII1907'],['OII3726','OII3729'],['NeIII3869'],['NeIII3967'],['HI4340'],['HI4861'],['OIII4959'],['OIII5007'],['HI6563']]

# lines = []


line = ['HI6563']

label = '-'.join(line)


labels = {}
labels['HI1216'] = r'Ly\alpha'
labels['CIII1909-CIII1907'] = r'\regular{[CII]1907,1909}'
labels['OII3726-OII3729'] = r'\regular{[OII]3726,3729}'
labels['NeIII3869'] = r'\regular{[NeIII]3869}'
labels['NeIII3967'] = r'\regular{[NeIII]3967}'
labels['HI4340'] = r'\regular{H\gamma}'
labels['HI4861'] = r'\regular{H\beta}'
labels['OIII5007'] = r'\regular{[OIII]5007}'
labels['OIII4959'] = r'\regular{[OIII]4959}'
labels['HI6563'] = r'\regular{H\alpha}'




# NOTE: SFR is implicitly 1

log10tmax = 7.

for SPS_IMF, ls in zip(SPS_models,line_styles):

    grid = pickle.load(open('../../BuildGrid/Z/grids/'+SPS_IMF+'/lines.p','rb'), encoding='latin1')


    line_luminosities = []

    for iZ, Z in enumerate(grid['Z']):

    
        line_luminosity = 0.0

        totSF = 0.0

        for ia, log10age in enumerate(np.arange(6.,log10tmax+0.1,0.1)):

            # --- determine SF activity in bin

            SF = 10**log10age - 10**(log10age-0.1) 
            if log10age==6.0: SF = 10**6.0     
            totSF += SF

            for l in line:
        
                l = l.encode('UTF-8')
        
                line_luminosity += SF*10**grid[l]['luminosity'][ia, iZ]


        line_luminosities.append(line_luminosity)
        
        
        
    ax.plot(grid['log10Z'], np.log10(np.array(line_luminosities)), label = SPS_labels[SPS_IMF], c='0.5', ls=ls)


        
# ax.legend(loc = 'center left', fontsize=7, bbox_to_anchor = (1.0, 0.5))

ax.legend(loc = 'lower left', fontsize=6, handlelength=2.5, labelspacing=0.2)


# ax.set_xlim([1.,3.])
ax.set_ylim([41., 42.])


ax.set_xlabel(r'$\rm\log_{10}(Z)$')
ax.set_ylabel(r'$\rm\log_{10}(L/erg\ s^{-1})$')


ax.axhline(41.27, c='k', alpha=0.1, lw=3)
ax.text(-5, 41.285, r'Kennicutt & Evans (2010)', fontsize = 6, color='k')


figname = 'figures/SFR_'+label+'_Z_IMF.pdf'
print(figname)
fig.savefig(figname)

# ---- define plot



