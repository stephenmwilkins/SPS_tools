
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


create_HDF5 = False

if create_HDF5:

    # --------------------------------
    # --- first of all make a simple copy of the Python structure and save this

    # --- create HDF5 file
    hf = h5py.File(f'{data_dir}/{SPS}/{IMF}/stellar.h5', 'w')

    # --- simply mirror the python data structure
    for k in p.keys():
        hf[k] = p[k]

    # --- also create a version where each SED is saved separately
    for iZ, Z in enumerate(p['Z']):
        for ia, log10age in enumerate(p['log10age']):
            hf[f'{Z}/{round(log10age,1)}/L_nu'] = p['L_nu'][ia, iZ]
            hf[f'{Z}/{round(log10age,1)}/fraction_remaining'] = p['fraction_remaining'][ia, iZ]


    hf.flush()
    hf.close()

os.system(f'cp {data_dir}/{SPS}/{IMF}/stellar.h5 {data_dir}/{SPS}_{IMF}.h5')
