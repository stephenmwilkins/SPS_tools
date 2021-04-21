
# --- this is designed to convert my pickle files to HDF5




import h5py


# where the data is
data_dir = '/Users/stephenwilkins/Dropbox/Research/data/SPS/stellar/1.0/'

# the SPS model and version
SPS = 'BPASSv2.2.1.binary'

# the initial mass function (IMF)
IMF = 'Chabrier_300'




# --- open HDF5 file

hf = h5py.File(f'{data_dir}/{SPS}/{IMF}/stellar.h5', 'r')


# --- explore datasets

for dset in hf: print(dset)

hf.visit(lambda name: print(name))


# ---

# print(hf['fraction_remaining'][:])
# print(hf['fraction_remaining'][:].shape)
