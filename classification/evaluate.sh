#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/tf

# set environment variables
export constrained=$1
export channel=$1
export constrained=$2
export modelpath=$3
export model=$4

mkdir -p /scratch/pahwagne/score_trees/$channel
python /work/pahwagne/RDsTools/classification/temp_evaluate.py
echo "copying files ..." 
xrdcp /scratch/pahwagne/score_trees/$channel/* root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/
echo "done!" 
rm -r /scratch/pahwagne/score_trees/$channel
