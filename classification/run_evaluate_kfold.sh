#!/bin/bash

#given by command line
export constrained=$1
export modelpath=$2

#remove "test_" from string
modelpath=${modelpath:5}

if [[ $constrained=="True" ]]; then 
  flag="_cons"
else 
  flag="_unc"
fi

echo $constrained
echo $modelpath

#prepare target diretcories if non existent
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/data_$modelpath$flag/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_$modelpath$flag/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/hb_$modelpath$flag/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/b0_$modelpath$flag/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/bplus_$modelpath$flag/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/bs_$modelpath$flag/

sbatch -p short -o log_data.txt  -e err_data.txt --mem=5G     evaluate_kfold.sh data       $constrained $modelpath 1 
sbatch -p short -o log_sig.txt   -e err_sig.txt --mem=5G      evaluate_kfold.sh sig        $constrained $modelpath 100 #for hammer weighting 
sbatch -p short -o log_hb.txt    -e err_hb.txt                evaluate_kfold.sh hb         $constrained $modelpath 1 
sbatch -p short -o log_bs.txt    -e err_bs.txt                evaluate_kfold.sh bs         $constrained $modelpath 1 
sbatch -p short -o log_b0.txt    -e err_b0.txt                evaluate_kfold.sh b0         $constrained $modelpath 1 
sbatch -p short -o log_bplus.txt -e err_bplus.txt             evaluate_kfold.sh bplus      $constrained $modelpath 1 

#interactive running!
#./evaluate_kfold.sh data  $constrained $modelpath 100 #interactvie running 
#./evaluate_kfold.sh sig   $constrained $modelpath 100 #interactvie running
#./evaluate_kfold.sh hb    $constrained $modelpath 100 #interactvie running 
#./evaluate_kfold.sh bplus $constrained $modelpath 100 #interactvie running 
#./evaluate_kfold.sh b0    $constrained $modelpath 100 #interactvie running 
#./evaluate_kfold.sh bs    $constrained $modelpath 100 #interactvie running 
