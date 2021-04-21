
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import SPS_tools.cloudy.abundances as a

print(a.abundances(0.02))
