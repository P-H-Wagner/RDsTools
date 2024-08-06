#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/tf

# set environment variables
export channel=$1
echo $channel

mkdir -p /scratch/pahwagne/score_trees
python /work/pahwagne/RDsTools/classification/temp_evaluate.py
echo "copying files ..." 
xrdcp /scratch/pahwagne/score_trees/* root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/
echo "done!" 
rm -r /scratch/pahwagne/score_trees
