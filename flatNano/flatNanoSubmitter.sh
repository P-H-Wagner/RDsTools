#!/bin/bash
cd /work/pahwagne/RDsTools/flatNano 
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/hammer3p8
mkdir -p /scratch/pahwagne/flatNano/
python flatNanoProducer.py
xrdcp -r /scratch/pahwagne/flatNano root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/
rm -r /scratch/pahwagne/flatNano
