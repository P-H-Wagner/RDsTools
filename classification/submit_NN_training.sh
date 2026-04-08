#!/bin/bash

dt=$1

#get the nr of folds from this model by parsing the corresponding dt folder

pnfs="/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nn_training/$dt/"

max_fold=$(ls $pnfs/xx_train_*.pck | sed -E 's/.*_([0-9]+)\.pck/\1/' | sort -n | tail -1)
echo "====> This model has $max_fold folds"

for i in $(seq 0 $max_fold); do

#write sh file to submit the training
cat << EOF > nn_training/$dt/submit_fold_$i.sh
#!/bin/bash

dt=$dt
fold=$i

python train_on_single_fold.py --datetime=$dt --fold=$i
EOF

#make it executable
chmod +x nn_training/$dt/submit_fold_${i}.sh

#now submit it to slurm
sbatch --job-name TRAIN_$i -p standard -o nn_training/$dt/out_$i.txt -e nn_training/$dt/err_$i.txt -n 1 --cpus-per-task=4 --mem-per-cpu=4000M nn_training/$dt/submit_fold_$i.sh
echo "sbatch --job-name TRAIN_$i -p short -o nn_training/$dt/out_$i.txt -e nn_training/$dt/err_$i.txt nn_training/$dt/submit_fold_$i.sh"

done

 
