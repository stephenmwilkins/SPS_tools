import streamlit as st
import numpy as np

# Use the non-interactive Agg backend, which is recommended as a
# thread-safe backend.
# See https://matplotlib.org/3.3.2/faq/howto_faq.html#working-with-threads.
import matplotlib as mpl
mpl.use("agg")

import matplotlib.pyplot as plt



import h5py
import matplotlib.pyplot as plt


# the SPS model and version and IMF
SPS = 'BPASSv2.2.1.binary'
IMF = 'Chabrier_300'


# --- open HDF5 file

hf = h5py.File(f'../sample_data/{SPS}_{IMF}_lite.h5', 'r')

lam = hf['lam'][:]

iZ = 8
ia = 10


iZ = st.sidebar.slider("Metallicity index", min_value=0, max_value=len(hf['Z'])-1, value=5, step=1)
st.sidebar.text(f"metallicity = {hf['Z'][iZ]}")
ia = st.sidebar.slider("Age index", min_value=0, max_value=len(hf['log10age'])-1, value=5, step=1)
st.sidebar.text(f"log10(age/yr) = {hf['log10age'][ia]}")


SED = hf['L_nu'][ia, iZ]

fig, ax = plt.subplots()

x = np.array([-1,1])

ax.plot(np.log10(lam), np.log10(SED))

ax.set_xlim([2.5, 4.0])
ax.set_ylim([16, 22])

# plt.title("y=mx+c")
plt.xlabel(r'$\log_{10}(\lambda/\AA)$')
plt.ylabel(r'$\log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1})$')


st.pyplot(fig)
