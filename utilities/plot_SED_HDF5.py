
# --- this is designed to convert my pickle files to HDF5


import h5py
import matplotlib.pyplot as plt

# where the data is
data_dir = '/Users/stephenwilkins/Dropbox/Research/data/SPS/stellar/1.0/'

# the SPS model and version
SPS = 'BPASSv2.2.1.binary'

# the initial mass function (IMF)
IMF = 'Chabrier_300'


# --- open HDF5 file

hf = h5py.File(f'{data_dir}/{SPS}/{IMF}/stellar.h5', 'r')

lam = hf['lam'][:]

iZ = 8
ia = 10

print(f"metallicity = {hf['Z'][iZ]}  log10(age/yr)={hf['log10age'][ia]}")


SED1 = hf['L_nu'][ia, iZ]


plt.plot(lam, SED1, alpha=0.2, lw=4, c='k')

# --- alternative approach
SED2 = hf[f"{hf['Z'][iZ]}/{hf['log10age'][ia]}/L_nu"]

plt.plot(lam, SED2, alpha=1, lw=1, c='k')


plt.show()
