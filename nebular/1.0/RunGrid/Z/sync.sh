#!/bin/bash

# push code/data to apollo ignoring coutputs and output

rsync -av --exclude='coutputs' --exclude='cinputs' --exclude='output' . apollo.hpc.susx.ac.uk:/research/astro/highz/SED/0.5/nebular/RunGrid/Z/ --delete



