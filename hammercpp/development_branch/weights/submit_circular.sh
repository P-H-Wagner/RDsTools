#!/bin/bash

#loop over files

export channel=$1
export input=$2
export target=$3
export nmax=$4
export prod=$5

echo $channel
echo $input
echo $target
echo $nmax
echo $prod

datetime=$(date +"%d_%m_%Y_%H_%M_%S")

if [ -z "$prod" ]; then
  echo "gen level MC, no production year provided, using <gen> as folder on pnfs"
  prod="gen"  # or "24", or some other fallback
fi

mkdir -p /work/pahwagne/RDsTools/hammercpp/development_branch/weights/${channel}/${datetime}/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/${prod}/${channel}_${target}_${datetime}/ 

make

# get variables from helper file by executing its main function
config_json=$(python /work/pahwagne/RDsTools/help/helper.py)
# get bash array out of json for the two signals (24 and 25)

if [ "$prod" = "24" ]; then
  echo "Getting 2024 production"
  eval "$(echo "$config_json" | jq -r '.sig_cons_24 | @sh "sig_cons=(\(.[]))"')"
  eval "$(echo "$config_json" | jq -r '.hb_cons_24  | @sh "hb_cons=(\(.[]))"')"
  echo "signal is: ${sig_cons[0]} and hb is: ${hb_cons[0]}"
fi

if [ "$prod" = "25" ]; then
  echo "Getting 2025 production"
  eval "$(echo "$config_json" | jq -r '.sig_cons_25 | @sh "sig_cons=(\(.[]))"')"
  eval "$(echo "$config_json" | jq -r '.hb_cons_25  | @sh "hb_cons=(\(.[]))"')"
  echo "signal is: ${sig_cons[0]} and hb is: ${hb_cons[0]}"
fi

echo "Getting unfiltered gen-level MC files"
dsmu_gen=$(echo "$config_json" | jq -r '.dsmu_gen')
dsmu_isgw2_gen=$(echo "$config_json" | jq -r '.dsmu_isgw2_gen')
dsstarmu_gen=$(echo "$config_json" | jq -r '.dsstarmu_gen')
dstau_gen=$(echo "$config_json" | jq -r '.dstau_gen')
dsstartau_gen=$(echo "$config_json" | jq -r '.dsstartau_gen')

echo $dsmu_gen
echo $dsmu_isgw2_gen
echo $dsstarmu_gen
echo $dstau_gen
echo $dsstartau_gen

if [ "$channel" = "signal" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/${sig_cons[0]}" # pre NN
    echo "signal"
fi

if [ "$channel" = "hb" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/${hb_cons[0]}" # pre NN
    echo "hb"
fi

if [ "$channel" = "dsmu" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/${dsmu_gen}"
    echo "dsmu"
fi

if [ "$channel" = "dsmu_isgw2" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/${dsmu_isgw2_gen}"
    echo "dsmu_isgw2"
fi

if [ "$channel" = "dsstarmu" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/${dsstarmu_gen}"
    echo "dsstarmu"
fi

if [ "$channel" = "dstau" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/${dstau_gen}"
    echo "dstau"
fi

if [ "$channel" = "dsstartau" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/${dsstartau_gen}"
    echo "dsstartau"
fi

#write submitter file for each chunk

counter=0

for file in "$path"/*; do

    echo "Hammering file: $file"
    #echo $path
    
    if (( nmax >= 0 && counter > nmax - 1)); then
        echo "n > nmax, breaking the loop"
        break
    fi

    echo $counter
    #write submitter file for each chunk
    submitter="/work/pahwagne/RDsTools/hammercpp/development_branch/weights/${channel}/${datetime}/submitter_${counter}.sh"
    cat >"$submitter" <<EOF
#!/bin/bash

export channel=\$1 #channel dsmu, dsmu_isgw2
export input=\$2   #input scheme
export target=\$3  #target scheme
export file=\$4    #file to hammer
export counter=\$5 #counter 

datetime=$datetime

eval "\$(conda shell.bash hook)"

cd /work/pahwagne/RDsTools/hammercpp/development_branch/weights/
mkdir -p /scratch/pahwagne/hammer/\$datetime
./hammer_circular \$channel \$input \$target \$file \$counter \$datetime
xrdcp /scratch/pahwagne/hammer/\$datetime/\${channel}_\${target}_\${counter}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/\${prod}/\${channel}_\${target}_\${datetime}/
rm /scratch/pahwagne/hammer/\$datetime/\${channel}_\${target}_\${counter}.root


EOF

    if [ -e "$file" ]; then
        cmd="sbatch -o ${channel}/${datetime}/out_${counter}.txt -e ${channel}/${datetime}/err_${counter}.txt -p standard -t 420 $submitter $channel $input $target $file $counter"
        sbatch -o ${channel}/${datetime}/out_${counter}.txt -e ${channel}/${datetime}/err_${counter}.txt -p standard -t 360 $submitter $channel $input $target $file $counter
        echo $cmd
    else
        echo "No files found in $path"
        break
    fi

    ((++counter))
done
