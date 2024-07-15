import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('id')
parser.add_argument('nDir')   # nr of directories on the t2 (counting from 0), f.e. nDir = 2 means we have: 0000,0001,0002 
parser.add_argument('nFiles') # nr of files to chain
args = parser.parse_args()

queue = 'short' 
#time = 60
nevents = -1

inputfiles = []

def list_files_gfal(path):

  # run gfal-ls and capture the output                   
  # stdout=subprocess.PIPE allows to pipe and save the output in the result variable rather than only printing it to the console
  # text=True gives more human readable format of output 

  result = subprocess.run(['gfal-ls', path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
  
  # Split the output by lines to get a list of files and directories
  files = result.stdout.splitlines()
  
  return files

for i in range(1): #range(int(args.nDir)):

  print(f"Chaining files from part {i}")

  path       = f"root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat//store/user/pahwagne/2024Jun28/ParkingBPH1/crab_{args.date_time}/{args.id}/000{i}/"
  files      = list_files_gfal(path)

  inputfiles = [path + f for f in files] #inlcude the path in the list :-)

os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/{args.date_time}")
os.makedirs(f"{args.date_time}/logs")
os.makedirs(f"{args.date_time}/errs")

for i,j in enumerate(range(0, len(inputfiles), int(args.nFiles))):

  #slice into lists of length args.nFiles
  fin = inputfiles[j:j+int(args.nFiles)]
  fout = f"/scratch/pahwagne/chained_data_{args.date_time}_chunk_{i}.root"

  #template
  temp = open("temp_chainer.py", "rt")
  #file to write to
  cfg = open(f"{args.date_time}/chainer.py","wt")
  
  for line in temp:
    if "HOOK_DATE_TIME"   in line: line = line.replace("HOOK_DATE_TIME", args.date_time)
    if "HOOK_FILE_IN"   in line: line = line.replace("HOOK_FILE_IN",     str(fin))
    if "HOOK_FILE_OUT"   in line: line = line.replace("HOOK_FILE_OUT",   fout)
    if "HOOK_CHUNK"   in line: line = line.replace("HOOK_CHUNK",         str(i))
  
    cfg.write(line)
  
  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/flatNano',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/{args.date_time}',
         f'python {args.date_time}/cfg_chunk_{i}.py',
         f'xrdcp /scratch/pahwagne/{args.date_time}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/{args.date_time}/{fout}',
         f'rm /scratch/pahwagne/{args.date_time}/{fout}',
         '',
     ])

  with open(f"{args.date_time}/submitter_chunk_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o {args.date_time}/logs/chunk_{i}.log',
        f'-e {args.date_time}/errs/chunk_{i}.err',
        '--mem=1000M',
        f'--job-name=FLAT_{i}',
        f'--time=5',
        f'{args.date_time}/submitter_chunk_{i}.sh',
     ])

  print(command_sh_batch)
  #os.system(command_sh_batch)







