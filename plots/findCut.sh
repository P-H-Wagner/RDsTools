#!/bin/bash

# loop over cuts
# there are no floats in bash, use awk
# produce 300 plots from 0 to 0.3 -> steps of permill

for x in $(awk 'BEGIN {for (i=0; i<=300; i++) printf "%.3f ", i/1000}')
do 
  echo "produce plots for cut $x"
  echo "sbatch -p short --time=20 -o logs/${x}_log.txt -o errs/${x}_err.txt submitPlots.sh $x"
  sbatch -p short --time=20 -o logs/${x}_log.txt -e errs/${x}_err.txt submitPlots.sh "$x"
  #sleep for 2s such that all 300 plotting directories have different names
  sleep 2
done
