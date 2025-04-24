import os
from glob import glob
from datetime import datetime
import sys 
import argparse

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
#parser.add_argument('channel') # sig or hb or data
parser.add_argument("-d", "--dt",     required = True, help = "Specify datetime")
parser.add_argument("-s", "--sb",     required = True, help = "Specify if bdt trained on 'left' sideband or both 'double'")
parser.add_argument("-n", "--nFiles", help = "Specify #files to process")
args = parser.parse_args()

if args.sb not in ["left", "double"]:
  raise ValueError ("Error: Not a valid key for --sb, please use 'left' or 'double'")
else: 
  sb = ""
  if args.sb == "double": sb += args.sb

print(sb)

model = "bdt_model_" + sb
dt_string = args.dt
######################################

#800 jobs per user on short queue
nMaxJobs = 1000

#default
queue = 'short' 
time = 30

filesPerJob = 30

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
for dt in data_cons_24:
  directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dt}/' #data bParking 2018 part D
  inputfiles +=  filesFromFolder(directory)

naming = 'data'

#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/"+dt_string)
os.makedirs(dt_string+"/logs")
os.makedirs(dt_string+"/errs")


for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_evaluate_bdt.py", "rt")
  #file to write to
  cfg = open(dt_string+"/eval_chunk_{0}.py".format(i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/{1}_chunk_{2}.root".format(dt_string,naming,i)  

  import pdb
  for line in temp:
    if   "HOOK_FILE_IN"  in line: 
      for f in fin:
        cfg.write(line.replace("HOOK_FILE_IN", f))
        
    elif "HOOK_FILE_OUT"  in line: cfg.write(line.replace("HOOK_FILE_OUT" , fout))
    elif "HOOK_DATE_TIME" in line: 
      line = line.replace("HOOK_DATE_TIME", dt_string)
      line = line.replace("HOOK_MODEL", model)
      cfg.write(line)
    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/classification/'+dt_string ,
         'mkdir -p /scratch/pahwagne/'+dt_string,
         'ls /scratch/pahwagne/',
         'python eval_chunk_{1}.py'.format(dt_string,i),
         'xrdcp /scratch/pahwagne/{0}/{1}_chunk_{2}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         'rm /scratch/pahwagne/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         '',
     ])

  with open("{0}/submitter_chunk_{1}.sh".format(dt_string,i), "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        '-o {0}/logs/chunk_{1}.log'.format(dt_string,i),
        '-e {0}/errs/chunk_{1}.err'.format(dt_string,i),
        '--mem=2500M',
        '--job-name=BDT_{0}'.format(i),
        '--time={0}'.format(time),
        '{0}/submitter_chunk_{1}.sh'.format(dt_string,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







