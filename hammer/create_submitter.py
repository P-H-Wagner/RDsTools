import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import subprocess

#example task: 240802_091822:pahwagne_crab_20240802_111807
#example command line: python create_submitter.py -f 0 -i  240802_091822 -p 2 20240802_111807 data t2

parser = argparse.ArgumentParser()

parser.add_argument('date_time') # second number sequence of crab task name
parser.add_argument('channel') # sig, data, hb, bs , ...
#parser.add_argument('site')  # t3 or t2 (from crab)
#parser.add_argument('-f','--folder')  # nr of directory on the t2, f.e. nDir = 2 means: 0002 
#parser.add_argument('-i','--crab_id') # first number sequence of crab task name 
#parser.add_argument('-p','--part')    #1,2,3,4,5 - specifies the part of the dataset (BPH1, BPH2, ..., BPH5) 
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


def list_files_gfal(path):

  # run gfal-ls and capture the output                   
  # stdout=subprocess.PIPE allows to pipe and save the output in the result variable rather than only printing it to the console
  # text=True gives more human readable format of output 

  result = subprocess.run(['gfal-ls', path], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

  # Split the output by lines to get a list of files and directories
  files = result.stdout.splitlines()

  # remove first element which is log (check this with printing it)
  files.remove('log')


  return files

#load flats
path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/" 
directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/"
inputfiles = os.listdir(directory)


#inputfiles = inputfiles[0:1] #debug
#print(inputfiles)

if not os.path.exists(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{args.date_time}"):

  os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{args.date_time}")
  os.makedirs(f"hammer_{args.date_time}/logs")
  os.makedirs(f"hammer_{args.date_time}/errs")

#inputfiles = inputfiles[0:1]

for fin in inputfiles:


  #get chunk nr
  pattern = r'([^_]+)_flatChunk_(\d+)\.root'
  match = re.search(pattern, fin)
  if match: chunkNr = match.group(2)      

  print(chunkNr)

  #template
  temp = open("hammer_dsstar_tau_temp.py", "rt")
  #file to write to
  cfg = open(f"hammer_{args.date_time}/cfg_chunk_{chunkNr}.py","wt")
  #file to save things
  fout = f"{naming}_hammerChunk_{chunkNr}.root"  

  for line in temp:
    if "HOOK_FILE_IN" in line:       line = line.replace("HOOK_FILE_IN", fin)
    if "HOOK_DATE_TIME" in line:   line = line.replace("HOOK_DATE_TIME", args.date_time)
    if "HOOK_FILE_OUT" in line:    line = line.replace("HOOK_FILE_OUT", f"/scratch/pahwagne/{args.date_time}/{fout}")

    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/hammer',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/{args.date_time}',
         f'python hammer_{args.date_time}/cfg_chunk_{chunkNr}.py',
         f'xrdcp /scratch/pahwagne/{args.date_time}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{args.date_time}/{fout}',
         f'rm /scratch/pahwagne/{args.date_time}/{fout}',
         '',
     ])

  with open(f"hammer_{args.date_time}/submitter_chunk_{chunkNr}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o hammer_{args.date_time}/logs/chunk_{chunkNr}.log',
        f'-e hammer_{args.date_time}/errs/chunk_{chunkNr}.err',
        #f'--mem=2000M',
        f'--job-name=HAMMER_{chunkNr}_{naming}',
        #f'--time=20',
        f'hammer_{args.date_time}/submitter_chunk_{chunkNr}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







