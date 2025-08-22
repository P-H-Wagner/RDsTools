#!/bin/bash

#given by command line
export modelpath=$1
export prod=$2

python create_submitter_kfold.py -m $modelpath -c sig   -p $prod 
python create_submitter_kfold.py -m $modelpath -c hb    -p $prod 
python create_submitter_kfold.py -m $modelpath -c data  -p $prod 
python create_submitter_kfold.py -m $modelpath -c bs    -p $prod 
python create_submitter_kfold.py -m $modelpath -c b0    -p $prod 
python create_submitter_kfold.py -m $modelpath -c bplus -p $prod 

