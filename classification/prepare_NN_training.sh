#!/bin/bash

#datetime is the first argument after .sh iteself
dt=$1
echo $dt
#dt=$(date +"%Y_%m_%d_%H_%M_%S")

scratch="/scratch/pahwagne/$dt/"
pnfs="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/$dt/"
pnfs_w_redir="root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/$dt/"

mkdir -p $scratch
mkdir -p $pnfs
ls -ld $scratch   

python prepare_NN_training.py --debug --datetime=$dt
#python prepare_NN_training.py --datetime=$dt
xrdcp -r  $scratch/* $pnfs_w_redir

rm -r $scratch
