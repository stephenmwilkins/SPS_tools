

# --- this demonstrates how to build a grid of Cloudy input files based on my standard set of parameters. An input file is created for every metallicity/age pair.


import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import SPS_tools.cloudy.BuildInput as BuildInput


# --- get the default set of parameters
SPS_params, cloudy_params = BuildInput.default_params()

# --- print the default parameters
for parameter_set in [SPS_params, cloudy_params]:
    for p, param in parameter_set.items():
        print(p, param)


# --- set the input stellar grid file location (will depend on where you put it)
SPS_params['stellar_grid_file'] = f"/Users/stephenwilkins/Dropbox/Research/data/SPS/stellar/1.0/{SPS_params['SPS']}/{SPS_params['IMF']}/stellar.p"


cloudy_params['output_dir'] = 'grids/basic/'

# --- make directory structure for the output files
if not os.path.exists(cloudy_params['output_dir']):
    os.makedirs(cloudy_params['output_dir'])


# --- initialise CloudyInput class. This allows you to access the input stellar grid that we use to generate the grid points

CI = BuildInput.CloudyInput(**SPS_params)


i = 1
for ia in range(len(CI.stellar_grid['log10age'])):
    for iZ in range(len(CI.stellar_grid['Z'])):

        print(ia, iZ)

        cloudy_params['output_file'] = f'{ia}_{iZ}' # simple
        cloudy_params['ia'] = ia
        cloudy_params['iZ'] = iZ

        cinput = CI.build_input(**cloudy_params)

        open(f"{cloudy_params['output_dir']}/{i}.in",'w').writelines(cinput)

        # --- for use with array jobs we also want to save each input sequentially

        i += 1




# strategy = 'simple'
# data_dir = '/research/astro/flare/data/SPS'
# output_dir = f'{data_dir}/nebular/CloudyOutputs/{strategy}/{cloudy_version}/{SPS}/{IMF}'
# output_filename = f'{ia}_{iZ}' # simple
# output_filename = f'{ia}_{iZ}_{log10n_H}_{log10U_S_0}_{CO}_{d2m}' # full
# output_file = f'{output_dir}/{output_filename}'
