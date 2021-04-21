
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import SPS_tools.cloudy.BuildInput as BuildInput


# --- get the default set of parameters
params = BuildInput.default_params()

# --- print the default parameters
for p, param in params.items():
    print(p, param)

# --- set the input stellar grid file location (will depend on where you put it)
params['stellar_grid_file'] = f"/Users/stephenwilkins/Dropbox/Research/data/SPS/stellar/1.0/{params['SPS']}/{params['IMF']}/stellar.p"

cinput = BuildInput.build_input(**params)


# strategy = 'simple'
# data_dir = '/research/astro/flare/data/SPS'
# output_dir = f'{data_dir}/nebular/CloudyOutputs/{strategy}/{cloudy_version}/{SPS}/{IMF}'
# output_filename = f'{ia}_{iZ}' # simple
# output_filename = f'{ia}_{iZ}_{log10n_H}_{log10U_S_0}_{CO}_{d2m}' # full
# output_file = f'{output_dir}/{output_filename}'
