import os
from glob import glob
from datetime import datetime
import sys 
import argparse

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
#parser.add_argument('channel') # sig or hb or data
parser.add_argument("-c",  "--channel"   ,   required = True, help = "Specify channel")
parser.add_argument("-dt", "--dt_string" ,   required = True, help = "Specify model name (dt)")
parser.add_argument("-s",  "--sel"       ,   required = True, help = "selection")
parser.add_argument("-n",  "--nFiles"    ,                    help = "Specify #files to process")
args = parser.parse_args()

now       = datetime.now()
#dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")
dt_string = args.dt_string
sel = args.sel
channel = args.channel

######################################

#800 jobs per user on short queue
#nMaxJobs = 1000

#default
queue = 'short' 
time = 15

filesPerJob = 1

######################################

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

if args.channel== "sig":
  directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_hammer_flatNano}/' #
  inputfiles +=  filesFromFolder(directory)
  naming = "sig"

if args.channel== "hb":
  directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{hb_flatNanos[0]}_fullBPark/' #
  inputfiles +=  filesFromFolder(directory)
  naming = "hb"

if args.channel== "data":
  #directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/' #

  for f in data_flatNanos:
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{f}_fullBPark/"
    print(f"adding file {base}")
    inputfiles +=  filesFromFolder(base) 

  #inputfiles +=  filesFromFolder(directory)

  naming = "data"

if "bph" in channel:

  #split by part
  for f in split_by_parts[channel]:
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.sel}/{f}_fullBPark/"
    print(f"====> adding file {base} ")
    inputfiles +=  filesFromFolder(base)

  naming = channel


print(f"total files {inputfiles}")
#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs( f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{naming}_{dt_string}/", exist_ok=True)
os.makedirs( f"{naming}_{dt_string}/logs" , exist_ok=True)
os.makedirs( f"{naming}_{dt_string}/errs" , exist_ok=True)


for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_evaluate_kfold_bdt.py", "rt")
  #file to write to
  cfg = open(f"{naming}_{dt_string}/eval_chunk_{i}.py","wt")
  #file to save things
  fout = f"/scratch/pahwagne/{dt_string}/{naming}_chunk_{i}.root"

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
         '#!/bin/bash',
         f"cd /work/pahwagne/RDsTools/classification/{naming}_{dt_string}/" ,
         f"mkdir -p /scratch/pahwagne/{dt_string}/",
         f"rm /scratch/pahwagne/{dt_string}/*" #clean folder
         f"ls /scratch/pahwagne/",
         f"python eval_chunk_{i}.py",
         f"xrdcp /scratch/pahwagne/{dt_string}/{naming}_chunk_{i}.root root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{naming}_{dt_string}/{naming}_chunk_{i}.root",
         f"rm /scratch/pahwagne/{dt_string}/{naming}_chunk_{i}.root",
         "",
     ])

  with open(f"{naming}_{dt_string}/submitter_chunk_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        f'--account=t3',
        f'-o {naming}_{dt_string}/logs/chunk_{i}.log',
        f'-e {naming}_{dt_string}/errs/chunk_{i}.err',
        f'--mem=6000M',
        f'--job-name=BDT_{i}',
        f'--time={time}',
        f'{naming}_{dt_string}/submitter_chunk_{i}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







