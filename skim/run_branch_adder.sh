#!/bin/bash

#export channel=$1

eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/hammer3p8
mkdir -p /scratch/pahwagne/adder/
cd /work/pahwagne/RDsTools/skim/
root -q -l 'branch_adder.cc("sig")'
xrdcp -r /scratch/pahwagne/adder/* root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/added/
rm -r /scratch/pahwagne

