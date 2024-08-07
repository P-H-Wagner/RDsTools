import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin

parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('channel') # which sig
parser.add_argument('site')  # t3 or t2 (from crab)
args = parser.parse_args()

queue = 'short' 
#time = 60
nevents = -1

#naming
if args.channel == 'sig': naming = 'all_signals'
elif args.channel == 'hb': naming = 'hb_inclusive'
elif args.channel == 'b0': naming = 'b0'
elif args.channel == 'bs': naming = 'bs'
elif args.channel == 'lambdab': naming = 'lambdab'
elif args.channel == 'bplus': naming = 'bplus'
else: naming = 'data'

#load nanos
if args.site == 't3':
  path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/" 
  directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/{args.date_time}/"
  inputfiles = os.listdir(directory)

else:
  path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/"
  directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/{args.date_time}/"
  inputfiles = os.listdir(directory)

#inputfiles = inputfiles[3307:3308]

os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}")
os.makedirs(f"flat_{args.date_time}/logs")
os.makedirs(f"flat_{args.date_time}/errs")

#inputfiles = inputfiles[0:5]

for fin in inputfiles:

  #get chunk nr
  pattern = r'([^_]+)_chunk_(\d+)\.root'
  match = re.search(pattern, fin)
  if match: chunkNr = match.group(2)

  print(chunkNr)

  #template
  temp = open("temp_flatNanoProducer.py", "rt")
  #file to write to
  cfg = open(f"flat_{args.date_time}/cfg_chunk_{chunkNr}.py","wt")
  #file to save things
  fout = f"{naming}_flatChunk_{chunkNr}.root"  

  for line in temp:
    if "HOOK_DATE_TIME" in line: line = line.replace("HOOK_DATE_TIME", args.date_time)
    if "HOOK_FILE_IN" in line:   line = line.replace("HOOK_FILE_IN", fin)
    if "HOOK_FILE_OUT" in line:  line = line.replace("HOOK_FILE_OUT", fout)
    if "HOOK_CHUNK" in line:     line = line.replace("HOOK_CHUNK", chunkNr)
    if "HOOK_CHANNEL" in line:   line = line.replace("HOOK_CHANNEL", naming)
    if "HOOK_PATH" in line:      line = line.replace("HOOK_PATH", path + '/' + args.date_time + '/' + fin)

    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/flatNano',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'cp nanoframe.py flat_{args.date_time}/nanoframe.py',
         f'mkdir -p /scratch/pahwagne/{args.date_time}',
         f'python flat_{args.date_time}/cfg_chunk_{chunkNr}.py',
         f'xrdcp /scratch/pahwagne/{args.date_time}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/{fout}',
         f'rm /scratch/pahwagne/{args.date_time}/{fout}',
         '',
     ])

  with open(f"flat_{args.date_time}/submitter_chunk_{chunkNr}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o flat_{args.date_time}/logs/chunk_{chunkNr}.log',
        f'-e flat_{args.date_time}/errs/chunk_{chunkNr}.err',
        f'--mem=2000M',
        f'--job-name=FLAT_{chunkNr}_{naming}',
        f'--time=30',
        f'flat_{args.date_time}/submitter_chunk_{chunkNr}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







