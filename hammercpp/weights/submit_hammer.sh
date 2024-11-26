#!/bin/bash

#loop over files

export channel=$1
export input=$2
export target=$3
export nmax=$4

mkdir -p /work/pahwagne/RDsTools/hammercpp/weights/${channel}/
mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/${channel}_${target}/ 
make

if [ "$channel" = "signal" ]; then
    path="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/26_07_2024_14_46_03"
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

#write submitter file for each chunk

counter=0


for file in "$path"/*; do

    
    if (( nmax >= 0 && counter > nmax - 1)); then
        echo "n > nmax, breaking the loop"
        break
    fi

    echo $counter
    #write submitter file for each chunk
    submitter="/work/pahwagne/RDsTools/hammercpp/weights/${channel}/submitter_${counter}.sh"
    cat >"$submitter" <<'EOF'
#!/bin/bash

export channel=$1 #channel dsmu, dsmu_isgw2
export input=$2 #input scheme
export target=$3 #target scheme
export file=$4 #file to hammer
export counter=$5 #counter 

eval "$(conda shell.bash hook)"

cd /work/pahwagne/RDsTools/hammercpp/weights
mkdir -p /scratch/pahwagne/hammer/
./hammer_temp $channel $input $target $file $counter 
xrdcp /scratch/pahwagne/hammer/${channel}_${target}_${counter}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/${channel}_${target}/
rm /scratch/pahwagne/hammer/${channel}_${target}_${counter}.root 

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
