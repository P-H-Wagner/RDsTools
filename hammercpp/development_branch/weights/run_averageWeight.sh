#!/bin/bash


datetime=$(date +"%d_%m_%Y_%H_%M_%S")

mkdir -p /work/pahwagne/RDsTools/hammercpp/development_branch/weights/${datetime}/

make

# get variables from helper file by executing its main function
config_json=$(python /work/pahwagne/RDsTools/help/helper.py)
# get bash array out of json for the two signals (24 and 25)

echo "Getting unfiltered gen-level MC files"
dsmu_to_bcl=$(echo "$config_json" | jq -r '.dsmu_to_bcl')
dstau_to_bcl=$(echo "$config_json" | jq -r '.dstau_to_bcl')
dsstarmu_to_bgl=$(echo "$config_json" | jq -r '.dsstarmu_to_bgl')
dsstartau_to_bgl=$(echo "$config_json" | jq -r '.dsstartau_to_bgl')

echo $dsmu_to_bcl
echo $dsstarmu_to_bg
echo $dstau_to_bcl
echo $dsstartau_to_bgl

cmake .
make
./getAverageWeight $dsmu_to_bcl  $dstau_to_bcl $dsstarmu_to_bgl $dsstartau_to_bgl $datetime
