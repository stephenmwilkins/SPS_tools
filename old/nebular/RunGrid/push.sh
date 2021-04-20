#!/bin/bash

# push code/data to apollo ignoring coutputs and output

rsync -av --exclude='coutputs' --exclude='cinputs' --exclude='output' . sw376@apollo.hpc.susx.ac.uk://research/astro/flare/utilities/SPS_tools/nebular/1.0/RunGrid/Z_refQ_wdust/ --delete
