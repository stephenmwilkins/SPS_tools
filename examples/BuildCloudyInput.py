

# --- this demonstrates how to build a Cloudy input file based on my standard set of parameters. There is also the option to run cloudy if you give it the location of the cloudy executable.


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


# --- now build the input file. This automatically saves it using the output_file parameter. In this case this is "test/test".
cinput = BuildInput.CloudyInput(**SPS_params).build_input(**cloudy_params)


# --- automatically run cloudy.
run_cloudy = False
if run_cloudy:
    os.system(f"/Users/stephenwilkins/Dropbox/Research/software/cloudy/c17.01/source/cloudy.exe < {cloudy_params['output_dir']}/{cloudy_params['output_file']}.in")
