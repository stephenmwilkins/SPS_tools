
import numpy as np
import pickle
import os

# --- Build stellar SED assuming continuous star formation activity.



SPS = 'P2'
IMFs = ['ModSalpeter_100']

SPS = 'BPASSv2.1.binary'
IMFs = ['ModSalpeter_100', 'ModSalpeter_300']

# SPS = 'BPASSv2.2.binary'
# IMFs = ['ModSalpeter_100', 'ModSalpeter_300']

# SPS = 'BPASSv2.2.1.binary'
# IMFs = ['1p0_100','1p0_300','ModSalpeter_100','ModSalpeter_300','1p7_100','1p7_300','Salpeter_100','Chabrier_100','Chabrier_300']


for IMF in IMFs:

    print(IMF)

    data_dir = '/Volumes/DATA/UtilitiesDATA/PopSynthesis/reformat/inputs/'+SPS

    log10ages = np.arange(6.0, 10.0, 0.1)

    lam = np.load(data_dir+'/lam.npy')

    s = [lam<20000.]

    nu = 3E8/(lam*1E-10)


    if SPS == 'P2': Zs = np.array([0.0004, 0.004, 0.008, 0.02, 0.05])  
    if SPS in ['BPASSv2.1.binary','BPASSv2.1.single','BPASSv2.2.binary','BPASSv2.2.single','BPASSv2.2.1.binary']: Zs = np.array([0.00001, 0.0001, 0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04])
    if SPS in ['BPASSv2.binary','BPASSv2.single']: Zs = np.array([0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04])

    log10Zs = np.log10(np.array(Zs))


    grid = {}

    grid['log10age'] = log10ages

    grid['log10Z'] = log10Zs
    grid['Z'] = np.array(Zs)

    grid['lam'] = lam[s]
    grid['L_nu'] = np.zeros( (len(log10ages), len(log10Zs), len(lam[s])) )
    grid['fraction_remaining'] = np.zeros( (len(log10ages), len(log10Zs)) )


    for iZ, Z in enumerate(Zs):
        
        mass_remaining = np.load(data_dir+'/'+IMF+'/'+'mass_'+str(Z)+'.npy')
        
        for iage, log10age in enumerate(log10ages):
          
            Lnu = np.load(data_dir+'/'+IMF+'/'+str(Z)+'_'+'{:.1f}'.format(log10age)+'.npy') # W Hz^-1 M_sol^-1
    
            Lnu *= 1E7 # erg Hz^-1 (M_sol yr^-1)^-1 
    
            grid['L_nu'][iage, iZ] = Lnu[s]
    
            grid['fraction_remaining'][iage, iZ] = mass_remaining[iage]
   
   
  


    outdir = '/'.join([SPS, IMF])

    if not os.path.exists('grids/'+SPS): os.mkdir('grids/'+SPS)   # --- make top model DIR

    if not os.path.exists('grids/'+outdir): os.mkdir('grids/'+outdir)   # --- make IMF model DIR

    pickle.dump(grid, open('grids/'+'/'+outdir+'/stellar.p','wb'))


