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
parser.add_argument('gen')  # hammer gen production for average weight? 
args = parser.parse_args()
 
queue = 'standard' 
time = 60
nevents = -1

#for now only hammer for signals
if args.channel == 'sig': 
  fname = sig_cons_pastNN 

#datetime of gen productions are specified in helper.py
if args.channel == "dstau":     gen_input = dstau_gen
if args.channel == "dsstartau": gen_input = dsstartau_gen

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
  path       = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/flat_{gen_input}/" 
  inputfiles = [path + f for f in os.listdir(path)]
  fname = args.channel

else:

  #load flats
  path       = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{fname}/" 
  inputfiles = [path + f for f in os.listdir(path)]
#print(inputfiles)

os.makedirs(f"hammer_{fname}/logs")
os.makedirs(f"hammer_{fname}/errs")

inputfiles = inputfiles[0:1]
#print(inputfiles)

for fin in inputfiles:


  #only simulated samples are divided into chunks
  #get chunk nr
  pattern = r'([^_]+)_flatChunk_(\d+)\.root'
  match = re.search(pattern, fin)
  if match: chunkNr = match.group(2)
  else: chunkNr = ''

  if args.gen == "True":
    template_file = "hammer_all_gen_temp.py"

  else:
    template_file = "hammer_all_temp.py"

  cfg_file      = f"hammer_{fname}/cfg_chunk_{chunkNr}.py"
  output_file   = f"{fname}_flatChunk_{chunkNr}.root"
  output_dir    = fname
  submit_file   = f"hammer_{fname}/submitter_chunk_{chunkNr}.sh"

  temp    = open(template_file,"rt")
  cfg     = open(cfg_file,"wt")
  fout    = output_file 
  foutdir = output_dir 

  for line in temp:
    if "HOOK_FILE_IN"    in line:       line = line.replace("HOOK_FILE_IN", fin)
    if "HOOK_MAX_EVENTS" in line:       line = line.replace("HOOK_MAX_EVENTS", str(nevents))
    if "HOOK_FILE_OUT"   in line:       line = line.replace("HOOK_FILE_OUT", f"/scratch/pahwagne/{foutdir}/{fout}")


    #print(f"/scratch/pahwagne/{fname}/{fout}")
    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/hammer',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/{foutdir}',
         f'python {cfg_file}',
         f'xrdcp /scratch/pahwagne/{foutdir}/{fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/{foutdir}/{fout}',
         f'rm /scratch/pahwagne/{foutdir}/{fout}',
         '',
     ])

  with open(submit_file, "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o hammer_{fname}/logs/chunk_{chunkNr}.log',
        f'-e hammer_{fname}/errs/chunk_{chunkNr}.err',
        #f'--mem=2000M',
        f'--job-name=HAMMER_{fname}',
        #f'--time={time}',
        f'{submit_file}',
        '--ntasks=1',
        '--cpus-per-task=8' 
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







