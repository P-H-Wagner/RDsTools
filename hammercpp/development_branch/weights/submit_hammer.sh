#!/bin/bash

#loop over files

export channel=$1
export input=$2
export target=$3
export nmax=$4
export prod=$5

datetime=$(date +"%d_%m_%Y_%H_%M_%S")

mkdir -p /work/pahwagne/RDsTools/hammercpp/development_branch/weights/${channel}/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/${channel}_${target}_${datetime}/ 

make

# get variables from helper file by executing its main function
config_json=$(python /work/pahwagne/RDsTools/help/helper.py)
# get bash array out of json for the two signals (24 and 25)

if [ "$prod" = "24" ]; then
  echo "Getting 2024 production"
  eval "$(echo "$config_json" | jq -r '.sig_cons_24 | @sh "sig_cons=(\(.[]))"')"
  eval "$(echo "$config_json" | jq -r '.hb_cons_24  | @sh "hb_cons=(\(.[]))"')"
fi

if [ "$prod" = "25" ]; then
  echo "Getting 2025 production"
  eval "$(echo "$config_json" | jq -r '.sig_cons_25 | @sh "sig_cons=(\(.[]))"')"
  eval "$(echo "$config_json" | jq -r '.hb_cons_25  | @sh "hb_cons=(\(.[]))"')"
fi

if [ "$channel" = "signal" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/${sig_cons[0]}" # pre NN
    echo "signal"
fi

if [ "$channel" = "hb" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/${hb_cons[0]}" # pre NN
    echo "signal"
fi

if [ "$channel" = "dsmu" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58"
    echo "dsmu"
fi

if [ "$channel" = "dsmu_isgw2" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/17_11_2024_17_26_37"
    echo "dsmu_isgw2"
fi

if [ "$channel" = "dsstarmu" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_09_43_21"
    echo "dsstarmu"
fi

if [ "$channel" = "dstau" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_13_50_26"
    echo "dstau"
fi

if [ "$channel" = "dsstartau" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/16_11_2024_09_45_34"
    echo "dsstartau"
fi

#write submitter file for each chunk

counter=0

for file in "$path"/*; do

    echo $file
    echo $path
    
    if (( nmax >= 0 && counter > nmax - 1)); then
        echo "n > nmax, breaking the loop"
        break
    fi

    echo $counter
    #write submitter file for each chunk
    submitter="/work/pahwagne/RDsTools/hammercpp/development_branch/weights/${channel}/submitter_${counter}.sh"
    cat >"$submitter" <<EOF
#!/bin/bash

export channel=\$1 #channel dsmu, dsmu_isgw2
export input=\$2   #input scheme
export target=\$3  #target scheme
export file=\$4    #file to hammer
export counter=\$5 #counter 

datetime=$datetime

eval "\$(conda shell.bash hook)"

cd /work/pahwagne/RDsTools/hammercpp/development_branch/weights
mkdir -p /scratch/pahwagne/hammer/
./hammer_temp \$channel \$input \$target \$file \$counter
xrdcp /scratch/pahwagne/hammer/\${channel}_\${target}_\${counter}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/\${prod}/\${channel}_\${target}_\${datetime}/
rm /scratch/pahwagne/hammer/\${channel}_\${target}_\${counter}.root


EOF

    if [ -e "$file" ]; then
        cmd="sbatch -o ${channel}/out_${counter}.txt -e ${channel}/err_${counter}.txt -p standard -t 420 $submitter $channel $input $target $file $counter"
        sbatch -o ${channel}/out_${counter}.txt -e ${channel}/err_${counter}.txt -p standard -t 360 $submitter $channel $input $target $file $counter
        echo $cmd
    else
        echo "No files found in $path"
        break
    fi

    ((++counter))
done
