#!/bin/bash

version=3.0

rsync -rtvu sw376@apollo.hpc.susx.ac.uk://research/astro/flare/data/SPS/nebular/$version/ /Users/stephenwilkins/research/data/SPS/nebular/$version/
