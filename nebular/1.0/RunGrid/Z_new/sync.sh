#!/bin/bash

# push code/data to apollo ignoring coutputs and output

rsync -av --exclude='coutputs' --exclude='cinputs' --exclude='output' . sw376@apollo.hpc.susx.ac.uk:/research/astro/highz/SED/0.5/nebular/RunGrid/Z_new/ --delete



