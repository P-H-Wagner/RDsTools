import os
from glob import glob
from datetime import datetime
import sys 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('channel')
args = parser.parse_args()

queue = 'short' 
time = 60
nevents = -1

#load nanos
directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/{args.date_time}/"

#get each file
inputfiles = os.listdir(directory)
#inputfiles = [directory + filename for filename in filenames ]

#inputfiles = inputfiles[0:1]

os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}")
os.makedirs(f"{args.date_time}/logs")
os.makedirs(f"{args.date_time}/errs")

for fin in inputfiles:

  chunkNr = fin[-7:-5]
  if '_' in chunkNr: chunkNr = chunkNr[-1]

  print(chunkNr)
  #template
  temp = open("temp_flatNanoProducer.py", "rt")
  #file to write to
  cfg = open(f"{args.date_time}/cfg_chunk_{chunkNr}.py","wt")
  #file to save things
  fout = f"{args.channel}_flatChunk_{chunkNr}.root"  

  for line in temp:
    if   "HOOK_DATE_TIME" in line: cfg.write(line.replace("HOOK_DATE_TIME", args.date_time))
    elif "HOOK_FILE_IN" in line: cfg.write(line.replace("HOOK_FILE_IN", fin))
    elif "HOOK_FILE_OUT" in line: cfg.write(line.replace("HOOK_FILE_OUT", fout))
    elif "HOOK_CHUNK" in line: cfg.write(line.replace("HOOK_CHUNK", chunkNr))
    elif "HOOK_CHUNK" in line: cfg.write(line.replace("HOOK_CHANNEL", args.channel))

    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/flatNano',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/flatNano/{args.date_time}',
         'python {args.date_time}/{cfg}',
         f'xrdcp /scratch/pahwagne/flatNano/{args.date_time}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/flatNano/{args.date_time}/{fout}',
         f'rm /scratch/pahwagne/flatNano/{args.date_time}/{fout}',
         '',
     ])

  with open(f"{args.date_time}/submitter_chunk_{chunkNr}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o {args.date_time}/logs/chunk_{chunkNr}.log',
        f'-e {args.date_time}/errs/chunk_{chunkNr}.err',
        '--mem=2500M',
        f'--job-name=FLAT_{chunkNr}',
        f'--time={time}',
        f'{args.date_time}/submitter_chunk_{chunkNr}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







