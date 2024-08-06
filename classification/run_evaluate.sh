#!/bin/bash

sbatch -p short -o log_data.txt  -e err_data.txt  --mem=5000M evaluate.sh data
sbatch -p short -o log_sig.txt   -e err_sig.txt               evaluate.sh sig
sbatch -p short -o log_hb.txt    -e err_hb.txt                evaluate.sh hb
sbatch -p short -o log_bs.txt    -e err_bs.txt                evaluate.sh bs
sbatch -p short -o log_b0.txt    -e err_b0.txt                evaluate.sh b0
sbatch -p short -o log_bplus.txt -e err_bplus.txt             evaluate.sh bplus

