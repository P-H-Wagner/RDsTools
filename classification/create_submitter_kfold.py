import os
from glob import glob
from datetime import datetime
import sys 
import argparse

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
#parser.add_argument('channel') # sig or hb or data
parser.add_argument("-m", "--model",     required = True, help = "Specify NN model in format 'test_28May2025_16h33m29s'")
parser.add_argument("-c", "--channel",   required = True, help = "Specify which channel to evalue (sig, hb, ..)")
parser.add_argument("-p", "--prod",   required = True, help = "Specify production year '24' or '25'")
parser.add_argument("-n", "--nFiles", help = "Specify #files to process")
args = parser.parse_args()

#remove test_ from model
modelpath = args.model[5:]
channel   = args.channel
prod      = args.prod

######################################

#800 jobs per user on short queue
#nMaxJobs = 1000

#default
queue = 'short' 
time = 30

filesPerJob = 3

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

# Input files

inputfiles = []
if prod == "24":

  baseline_selection = base_wout_tv_24
 
  if args.channel == "data":
    #bdt is evaluated on the skimmed datasets :) 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_24}/"
    inputfiles +=  filesFromFolder(base) 
 
  if args.channel == "sig": 
    #hammer is evaluated on the skimmed datasets :) 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{sig_cons_hammer_24}/" 
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "hb": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_24[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "bs": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_24[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "b0": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_24[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "bplus": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_24[0]}/"
    inputfiles +=  filesFromFolder(base) 

elif prod == "25":

  baseline_selection = base_wout_tv_25
  # for now we only consider mu7
  baseline_selection += " && (mu7_ip4) "

  if args.channel == "data":
    #bdt is evaluated on the skimmed datasets :) 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/"
    inputfiles +=  filesFromFolder(base) 

  if args.channel == "sig": 
    #hammer is evaluated on the skimmed datasets :) 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{sig_cons_hammer_25}/" 
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "hb": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{hb_cons_25[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "bs": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bs_cons_25[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "b0": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{b0_cons_25[0]}/"
    inputfiles +=  filesFromFolder(base) 
  
  if args.channel == "bplus": 
    base = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{bplus_cons_25[0]}/"
    inputfiles +=  filesFromFolder(base) 
 
#####################################


if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{channel}_{modelpath}/", exist_ok=True)
os.makedirs(f"evaluation/{modelpath}/{args.channel}/"     , exist_ok=True)
os.makedirs(f"evaluation/{modelpath}/{args.channel}/logs" , exist_ok=True)
os.makedirs(f"evaluation/{modelpath}/{args.channel}/errs" , exist_ok=True)


for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_evaluate_kfold.py", "rt")
  #file to write to
  cfg = open(f"evaluation/{modelpath}/{args.channel}/eval_chunk_{i}.py", "wt")
  #file to save things
  fout = f"/scratch/pahwagne/score_trees/{args.channel}/{args.channel}_{modelpath}_flatChunk_{i}.root"  


  file_list = ', '.join(f'"{f}"' for f in fin)
  placeholders = {
    "HOOK_FILES":      file_list,
    "HOOK_MODELPATH":  modelpath,
    "HOOK_SELEC":      baseline_selection,
    "HOOK_CHANNEL":    args.channel,
    "HOOK_FILE_OUT":   fout,
    "HOOK_CHUNK":      str(i) 
  }

  for line in temp:
    for key, val in placeholders.items():
      line = line.replace(key, val)
    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         f'cd /work/pahwagne/RDsTools/classification/evaluation/{modelpath}/{args.channel}' ,
         f'mkdir -p /scratch/pahwagne/score_trees/{args.channel}/',
         'ls /scratch/pahwagne/',
         f'python eval_chunk_{i}.py',
         f'xrdcp /scratch/pahwagne/score_trees/{args.channel}/{args.channel}_{modelpath}_flatChunk_{i}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{args.channel}_{modelpath}/{args.channel}_{modelpath}_flatChunk_{i}.root',
         f'rm /scratch/pahwagne/score_trees/{args.channel}/{args.channel}_{modelpath}_flatChunk_{i}.root',
     ])

  with open(f"evaluation/{modelpath}/{args.channel}/submitter_chunk_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        f'-o evaluation/{modelpath}/{args.channel}/logs/chunk_{i}.log',
        f'-e evaluation/{modelpath}/{args.channel}/errs/chunk_{i}.err',
        '--mem=2000M',
        f'--job-name=EVAL_{i}',
        f'--time={time}',
        f'evaluation/{modelpath}/{args.channel}/submitter_chunk_{i}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







