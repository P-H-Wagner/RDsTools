#!/bin/bash

export var=$1
datetime=$(date +"%d_%m_%Y_%H_%M_%S")

mkdir -p /work/pahwagne/RDsTools/hammercpp/development_branch/weights/plots/${datetime}/

#run histModels.py to get plotting settings
#be sure the variable var is included in the models {} in histModels.py !
python /work/pahwagne/RDsTools/plots/histModels.py

# get variables from helper file by executing its main function
config_json=$(python /work/pahwagne/RDsTools/help/helper.py)
# get bash array out of json for the two signals (24 and 25)

echo "Getting unfiltered gen-level MC files"
sig_cons_hammer_25=$(echo "$config_json" | jq -r '.sig_cons_hammer_25')
averageWeightsYaml=$(echo "$config_json" | jq -r '.averageWeightsYaml')
base_wout_tv_25=$(echo "$config_json" | jq -r '.base_wout_tv_25')

# extend selection to mu7 and cos cut:
echo $base_wout_tv_25
base_wout_tv_25="${base_wout_tv_25} && (mu7_ip4) && (mu_is_global) && (ds_vtx_cosine_xy_pv > 0.8)"

echo $sig_cons_hammer_25
echo $averageWeightsYaml
echo $base_wout_tv_25

cmake .
make
./plotVariations $var $sig_cons_hammer_25 $averageWeightsYaml "${base_wout_tv_25}" $datetime
