#!/bin/bash

#loop over files

export channel=$1
export input=$2
export target=$3
export nmax=$4

datetime=$(date +"%d_%m_%Y_%H_%M_%S")

mkdir -p /work/pahwagne/RDsTools/hammercpp/weights/${channel}/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/${channel}_${target}_${datetime}/ 

make

if [ "$channel" = "signal" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_26Sep2024_07h46m21s_cons"
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
    submitter="/work/pahwagne/RDsTools/hammercpp/weights/${channel}/submitter_${counter}.sh"
    cat >"$submitter" <<EOF
#!/bin/bash

export channel=\$1 #channel dsmu, dsmu_isgw2
export input=\$2   #input scheme
export target=\$3  #target scheme
export file=\$4    #file to hammer
export counter=\$5 #counter 

datetime=$datetime

eval "\$(conda shell.bash hook)"

cd /work/pahwagne/RDsTools/hammercpp/weights
mkdir -p /scratch/pahwagne/hammer/
./checkHammer \$channel \$input \$target \$file \$counter
xrdcp /scratch/pahwagne/hammer/\${channel}_\${target}_\${counter}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/\${channel}_\${target}_\${datetime}/
rm /scratch/pahwagne/hammer/\${channel}_\${target}_\${counter}.root


EOF

    if [ -e "$file" ]; then
        cmd="sbatch -o ${channel}/out_${counter}.txt -e ${channel}/err_${counter}.txt -p short $submitter $channel $input $target $file $counter"
        sbatch -o ${channel}/out_${counter}.txt -e ${channel}/err_${counter}.txt -p short $submitter $channel $input $target $file $counter
        echo $cmd
    else
        echo "No files found in $path"
        break
    fi

    ((++counter))
done
