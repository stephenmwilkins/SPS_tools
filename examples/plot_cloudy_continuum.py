
import os
import sys

import numpy as np
import pickle
import yaml

import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import SPS_tools.cloudy.read_output as read_output

# --- defines the plot style. Can be safely commented out
import FLARE.plt as fplt



model = 'test/test'

# ---
output_dir, output_file = model.split('/')

# --- open yaml parameters file (contains SPS, IMF, ia, iZ, etc.)

with open(f'{model}.yaml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    params = yaml.load(file, Loader=yaml.FullLoader)







# --- explain why this is needed
Norm = params['log10Q_orig'] - params['log10Q'] # ARGG... no idea why + 1

# --- print parameters
for k, v in params.items():
    print(k, v)


# --- read direct cloudy spectra

lam, nu, incident, transmitted, nebular_continuum, total, linecont = read_output.continuum(output_dir, output_file)


for xlim, label in zip([False, [2.5, 4.0]], ['all', 'UVOptNIR']):

    fig = plt.figure(figsize = (4,4))

    left  = 0.15
    bottom = 0.15
    height = 0.8
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))


    for continuum_type in ['incident', 'transmitted', 'nebular_continuum', 'total']:

        ax.plot(np.log10(lam), np.log10(locals()[continuum_type]) + Norm, label = continuum_type, alpha = 0.5)

    stellar_grid = pickle.load(open(params['stellar_grid_file'],'rb'))

    ax.plot(np.log10(stellar_grid['lam']), np.log10(stellar_grid['L_nu'][int(params['ia']), int(params['iZ'])]), label = 'input incident', ls = '--', alpha = 0.5)

    if xlim:
        ax.set_xlim(xlim)
    ax.set_ylim([19,22])

    ax.set_xlabel(r'$\rm \log_{10}(\lambda/\AA m)$')
    ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\ M_{\odot}^{-1})$')

    ax.legend(loc = 'lower left', fontsize = 8)


    # info = '\n'.join([f'{k}={v}' for k, v in params.items()])
    # ax.text(0.025, 0.975, info, size=6, va="top", ha="left", multialignment="left", transform=ax.transAxes)



    fig.savefig(f'{model}_cont_{label}.pdf')

    fig.clf()
