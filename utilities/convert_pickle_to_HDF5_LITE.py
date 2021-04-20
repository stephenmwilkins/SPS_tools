
# --- this is designed to convert my pickle files to HDF5




import pickle
import h5py
import numpy as np
import os


# where the data is
data_dir = '/Users/stephenwilkins/Dropbox/Research/data/SPS/stellar/1.0/'

# the SPS model and version
SPS = 'BPASSv2.2.1.binary'

# the initial mass function (IMF)
IMF = 'Chabrier_300'


p = pickle.load(open(f'{data_dir}/{SPS}/{IMF}/stellar.p','rb'))

print(p.keys())

# print(list(map(str, p['Z'])))


p['log10age'] = np.round(p['log10age'], 1)


create_HDF5 = True

if create_HDF5:

    # --------------------------------
    # --- first of all make a simple copy of the Python structure and save this

    # --- create HDF5 file
    hf = h5py.File(f'{data_dir}/{SPS}_{IMF}_lite.h5', 'w')

    # --- simply mirror the python data structure
    for k in ['Z','log10age']:
        hf[k] = p[k]

    lam = p['lam']
    lam = lam[:10000]
    lam_rebinned = (lam[0::4] + lam[1::4] + lam[2::4] + lam[3::4])/4
    hf['lam'] = lam_rebinned

    hf['L_nu'] = np.ones((len(p['log10age']), len(p['Z']), 2500))

    # --- also create a version where each SED is saved separately
    for iZ, Z in enumerate(p['Z']):
        for ia, log10age in enumerate(p['log10age']):

            L_nu = p['L_nu'][ia, iZ][:10000]
            L_nu_rebinned = (L_nu[0::4] + L_nu[1::4] + L_nu[2::4] + L_nu[3::4])/4

            hf['L_nu'][ia, iZ, :] = L_nu_rebinned


    hf.flush()
    hf.close()
