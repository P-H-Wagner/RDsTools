#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/hammer3p8
python /work/pahwagne/RDsTools/comb/plotNtuple.py
