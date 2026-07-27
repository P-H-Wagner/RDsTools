import os
from glob import glob
from datetime import datetime
import sys 
import argparse

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nFiles", help = "Specify #files to process")
parser.add_argument("--part",        help = "specify part to process (i.e. 1,2,3, ...)")
args = parser.parse_args()

now       = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

dt_string = bdt_model #use the same name as the model!


######################################

#800 jobs per user on short queue
#nMaxJobs = 1000

#default
queue = 'short' 
time = 15 

filesPerJob = 1

######################################

#queue = 'standard'
#time = 60
nevents = -1
if args.nFiles: 
  nFiles = int(args.nFiles)
else: 
  nFiles = -1

######################################

def filesFromFolder(direc):
  filenames = os.listdir(direc)
  return [direc + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line[0:-1] for line in dataFiles] #-2 to remove the newline character \n
  return [filename for filename in filenames ]

#######################################

# Input source

#loop over all bph parts
inputfiles = []
if args.part:
  files_data = f"data_bph{args.part}_{nn_model}"
else:
  files_data = f"data_{nn_model}" 
#
print(f"====> evaluating BDT on {files_data}")
directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{files_data}/' #
inputfiles +=  filesFromFolder(directory)

#for f in split_by_parts[f"data_bph{args.part}"]:
#
#  directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/with_iso/{f}_fullBPark/' 
#  print(f"Adding files of {directory}")
#  inputfiles +=  filesFromFolder(directory)

naming = f'data_bph{args.part}'

#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

outdir = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{dt_string}/{naming}_{dt_string}/"
os.makedirs(outdir, exist_ok=True)
os.makedirs(f"{dt_string}/{naming}_{dt_string}/logs" , exist_ok=True)
os.makedirs(f"{dt_string}/{naming}_{dt_string}/errs" , exist_ok=True)


for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_evaluate_bdt.py", "rt")
  #file to write to
  cfg = open(f"{dt_string}/{naming}_{dt_string}/eval_chunk_{i}.py","wt")
  #file to save things
  fout = f"/scratch/pahwagne/{naming}_{dt_string}/{naming}_chunk_{i}.root"  

  import pdb
  for line in temp:
    if   "HOOK_FILE_IN"  in line: 
      for f in fin:
        cfg.write(line.replace("HOOK_FILE_IN", f))
        
    elif "HOOK_FILE_OUT"  in line: cfg.write(line.replace("HOOK_FILE_OUT" , fout)     )
    elif "HOOK_DATE_TIME" in line: cfg.write(line.replace("HOOK_DATE_TIME", dt_string))
    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         f'#!/bin/bash',
         f'eval "$(conda shell.bash hook)"',
         f'conda activate /work/pahwagne/environments/gpu_bdt',
         f'cd /work/pahwagne/RDsTools/classification/{dt_string}/{naming}_{dt_string}',
         f'mkdir -p /scratch/pahwagne/{naming}_{dt_string}/',
         f'ls /scratch/pahwagne/',
         f'python eval_chunk_{i}.py',
         f'xrdcp {fout}  root://t3dcachedb03.psi.ch:1094//{outdir}/{naming}_chunk_{i}.root',
         f'rm {fout}',
         f'',
     ])

  with open(f"{dt_string}/{naming}_{dt_string}/submitter_chunk_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)

  command_sh_batch = ' '.join([

        f'sbatch',
        f'-p '+queue,
        f'--account=t3',
        f'-o {dt_string}/{naming}_{dt_string}/logs/chunk_{i}.log',
        f'-e {dt_string}/{naming}_{dt_string}/errs/chunk_{i}.err',
        #f'--mem=6000M',
        f'--job-name=BDT_{i}',
        f'--time={time}',
        f'{dt_string}/{naming}_{dt_string}/submitter_chunk_{i}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







