import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import subprocess
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

#example task: 240802_091822:pahwagne_crab_20240802_111807
#example command line: python create_submitter.py -f 0 -i  240802_091822 -p 2 20240802_111807 data t2

parser = argparse.ArgumentParser()

parser.add_argument('channel') # sig, data, hb, bs , ...
parser.add_argument('-g','--gen')  # gen production 
parser.add_argument('-dt','--date_time') # second number sequence of crab task name
#parser.add_argument('-f','--folder')  # nr of directory on the t2, f.e. nDir = 2 means: 0002 
#parser.add_argument('-i','--crab_id') # first number sequence of crab task name 
#parser.add_argument('-p','--part')    #1,2,3,4,5 - specifies the part of the dataset (BPH1, BPH2, ..., BPH5) 
args = parser.parse_args()
 
#print(args.date_time)

queue = 'standard' 
time = 60
nevents = -1

#naming (for now only hammer for signals)
if args.channel == 'sig': 
  naming = 'sig' 
  fname = sig_cons_pastNN 


#get date_time of evaluated signal root tree


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

if args.gen == "True":

  #load flats
  path       = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/" 
  #directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{fname}"                     else:
  #inputfiles = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/{fname}.root"  #only one file!
  #print(path)

  inputfiles = [path + f for f in os.listdir(path)]
  fname = args.channel
  naming = args.channel

else:

  #load flats
  path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/" 
  #directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{fname}"
  inputfiles = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{fname}.root"  #only one file!

  inputfiles = [inputfiles]
#print(inputfiles)

os.makedirs(f"hammer_{fname}/logs")
os.makedirs(f"hammer_{fname}/errs")

#inputfiles = inputfiles[0:1]
#print(inputfiles)

for fin in inputfiles:

  #print(fin)
  #get chunk nr
  pattern = r'([^_]+)_flatChunk_(\d+)\.root'
  match = re.search(pattern, fin)
  if match: chunkNr = match.group(2)
  else: chunkNr = '0'

  if args.gen == "True":
    temp = open("hammer_all_gen_temp.py", "rt")

  else:
    temp = open("hammer_all_temp.py", "rt")
  #file to write to
  cfg = open(f"hammer_{fname}/cfg_chunk_{chunkNr}.py","wt")
  #file to save things
  fout = f"{naming}_flatChunk_{chunkNr}.root"

  for line in temp:
    if "HOOK_FILE_IN"    in line:       line = line.replace("HOOK_FILE_IN", fin)
    if "HOOK_MAX_EVENTS" in line:       line = line.replace("HOOK_MAX_EVENTS", str(nevents))
    if "HOOK_FILE_OUT"   in line:    line = line.replace("HOOK_FILE_OUT", f"/scratch/pahwagne/{fname}/{fout}")


    #print(f"/scratch/pahwagne/{fname}/{fout}")
    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/hammer',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/{fname}',
         f'python hammer_{fname}/cfg_chunk_{chunkNr}.py',
         f'xrdcp /scratch/pahwagne/{fname}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{fout}',
         f'rm /scratch/pahwagne/{fname}/{fout}',
         '',
     ])

  with open(f"hammer_{fname}/submitter_chunk_{chunkNr}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o hammer_{fname}/logs/chunk_{chunkNr}.log',
        f'-e hammer_{fname}/errs/chunk_{chunkNr}.err',
        #f'--mem=2000M',
        f'--job-name=HAMMER_{naming}',
        #f'--time={time}',
        f'hammer_{fname}/submitter_chunk_{chunkNr}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







