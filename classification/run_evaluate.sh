#!/bin/bash


export constrained=$1
export modelpath=$2
export model=$3

echo $constrained
echo $modelpath
echo $model


sbatch -p short -o log_data.txt  -e err_data.txt  --mem=5000M evaluate.sh  data   $constrained $modelpath $model
sbatch -p short -o log_sig.txt   -e err_sig.txt               evaluate.sh sig   $constrained $modelpath $model
sbatch -p short -o log_hb.txt    -e err_hb.txt                evaluate.sh hb    $constrained $modelpath $model  
sbatch -p short -o log_bs.txt    -e err_bs.txt                evaluate.sh bs    $constrained $modelpath $model
sbatch -p short -o log_b0.txt    -e err_b0.txt                evaluate.sh b0    $constrained $modelpath $model
sbatch -p short -o log_bplus.txt -e err_bplus.txt             evaluate.sh bplus $constrained $modelpath $model

