#!/bin/bash


datetime=$(date +"%d_%m_%Y_%H_%M_%S")

mkdir -p /work/pahwagne/RDsTools/hammercpp/tests/${datetime}/

make

# get variables from helper file by executing its main function
config_json=$(python /work/pahwagne/RDsTools/help/helper.py)
# get bash array out of json for the two signals (24 and 25)

echo "Getting unfiltered gen-level MC files"
dsmu_to_isgw2=$(echo "$config_json" | jq -r '.dsmu_to_isgw2')
dsmu_isgw2_to_cln=$(echo "$config_json" | jq -r '.dsmu_isgw2_to_cln')

echo $dsmu_to_isgw2
echo $dsmu_isgw2_to_cln

cmake .
make
./circularTest $dsmu_to_isgw2 $dsmu_isgw2_to_cln $datetime
