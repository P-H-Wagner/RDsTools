#!/bin/bash

dt=$(date +"%Y_%m_%d_%H_%M_%S")
echo $dt

here="/work/pahwagne/RDsTools/classification/nn_training/$dt/"
pnfs="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/$dt/"

mkdir -p $here
mkdir -p $pnfs

#queue submission
#sbatch -p short -o $here/out.log -e $here/err.err --mem=6000M prepare_NN_training.sh $dt 

#local
#prepare_NN_training.sh $dt 

