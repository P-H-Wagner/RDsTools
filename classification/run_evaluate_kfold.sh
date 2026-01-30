#!/bin/bash

#given by command line
export trigger=$1
export prod=$2

python create_submitter_kfold.py -t $trigger  -c sig   -p $prod  
python create_submitter_kfold.py -t $trigger  -c hb    -p $prod 
#python create_submitter_kfold.py -t $trigger  -c data  -p $prod 
#python create_submitter_kfold.py -m $modelpath -c bs    -p $prod 
#python create_submitter_kfold.py -m $modelpath -c b0    -p $prod 
#python create_submitter_kfold.py -m $modelpath -c bplus -p $prod 

