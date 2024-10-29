#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/tf

# given by command line 
export channel=$1
export constrained=$2
export modelpath=$3
export nchunks=$4

if [[ $constrained=="True" ]]; then 
  flag="_cons"
else 
  flag="_unc"
fi

target="${channel}_${modelpath}${flag}"

mkdir -p /scratch/pahwagne/score_trees/$channel
python /work/pahwagne/RDsTools/classification/temp_evaluate_kfold.py
echo "copying files ..." 
xrdcp /scratch/pahwagne/score_trees/$channel/* root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/$target/
echo "done!" 
rm -r /scratch/pahwagne/score_trees/$channel
