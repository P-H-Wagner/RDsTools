#!/bin/bash

dt=$1

#get the nr of folds from this model by parsing the corresponding dt folder
pnfs="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/bdt_training/$dt/"

#prepare a folder to save all the outputs
mkdir -p ./bdt_outputs/bdt_${dt}

max_fold=$(ls $pnfs/X_train_*.pck | sed -E 's/.*_([0-9]+)\.pck/\1/' | sort -n | tail -1)
echo "====> This model has $max_fold folds"

for i in $(seq 0 $max_fold); do

#write sh file to submit the training
cat << EOF > bdt_training/$dt/submit_fold_$i.sh
#!/bin/bash
source /work/pahwagne/miniconda3/etc/profile.d/conda.sh
conda activate cpu_bdt

dt=$dt
fold=$i

python train_on_single_fold_BDT.py --datetime=$dt --fold=$i
EOF

#make it executable
chmod +x bdt_training/$dt/submit_fold_${i}.sh

#now submit it to slurm
sbatch --job-name TRAIN_$i -p standard -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --cpus-per-task=8 --mem-per-cpu=4000M bdt_training/$dt/submit_fold_$i.sh
#sbatch --job-name TRAIN_$i -p gpu -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --cpus-per-task=8 --mem-per-cpu=4000M bdt_training/$dt/submit_fold_$i.sh
#sbatch --job-name TRAIN_$i -p short -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --cpus-per-task=4 --mem-per-cpu=4000M bdt_training/$dt/submit_fold_$i.sh
#sbatch --job-name TRAIN_$i -p gpu -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --cpus-per-task=8 --mem-per-cpu=6000M --time=240 --account=gpu_gres bdt_training/$dt/submit_fold_$i.sh
#sbatch --job-name TRAIN_$i -p gpu -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --mem=20000M --time=240 --account=gpu_gres bdt_training/$dt/submit_fold_$i.sh
#echo "sbatch --job-name TRAIN_$i -p short -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt bdt_training/$dt/submit_fold_$i.sh"

echo "sbatch --job-name TRAIN_$i -p gpu -o bdt_training/$dt/out_$i.txt -e bdt_training/$dt/err_$i.txt -n 1 --cpus-per-task=4 --mem-per-cpu=4000M --time=240 --account=gpu_gres bdt_training/$dt/submit_fold_$i.sh"
done

 
