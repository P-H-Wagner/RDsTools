#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /work/pahwagne/environments/hammer3p8

#cut is the cut on the discriminator, score5

export cut="$1"
echo "cut given to submitPlots.sh is: $cut"
python plotter.py --fit=cons --hb=old --nn=past --hammer=true --hammer_sys=false --prod=25 --bdt=true --debug --cut=$cut --findcut

