import os
from glob import glob
from datetime import datetime
import sys 
import argparse

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
#parser.add_argument('channel') # sig or hb or data
parser.add_argument("-p", "--prod",   required = True, help = "Specify production year '24' or '25'")
#parser.add_argument("-t", "--trigger",required = True, help = "Specify '7' or '9' to specify trigger menu")
parser.add_argument("-n", "--nFiles", help = "Specify #files to process")
args = parser.parse_args()

#if args.trigger not in ["7", "9"]:
#  raise ValueError ("Error: Not a valid key for --trigger, please use '7' or '9'")
#else: trig = args.trigger

folder = ""
now       = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

######################################

#800 jobs per user on short queue
#nMaxJobs = 1000

#default
queue = 'short' 
time = 30

filesPerJob = 6

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
if args.prod == "24":
  for dt in data_cons_24:
    directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dt}/' #
    inputfiles +=  filesFromFolder(directory)
else:
  for dt in data_cons_25:
    directory = f'/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dt}/' #
    inputfiles +=  filesFromFolder(directory)

naming = 'data'

#####################################

if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/"+ folder + "/" +dt_string, exist_ok=True)
os.makedirs(dt_string + "/" + folder  + "/logs" , exist_ok=True)
os.makedirs(dt_string + "/" + folder  + "/errs" , exist_ok=True)


for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_evaluate_bdt.py", "rt")
  #file to write to
  cfg = open(dt_string+"/{0}/eval_chunk_{1}.py".format(folder,i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/{1}_chunk_{2}.root".format(dt_string,naming,i)  

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
         'cd /work/pahwagne/RDsTools/classification/'+dt_string + "/" + folder ,
         'mkdir -p /scratch/pahwagne/'+dt_string,
         'ls /scratch/pahwagne/',
         'python eval_chunk_{1}.py'.format(dt_string,i),
         'xrdcp /scratch/pahwagne/{1}/{2}_chunk_{3}.root root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{0}/{1}/{2}_chunk_{3}.root'.format(folder,dt_string,naming,i),
         'rm /scratch/pahwagne/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         '',
     ])

  with open("{0}/{1}/submitter_chunk_{2}.sh".format(dt_string,folder,i), "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        '-o {0}/{1}/logs/chunk_{2}.log'.format(dt_string,folder,i),
        '-e {0}/{1}/errs/chunk_{2}.err'.format(dt_string,folder,i),
        '--mem=3000M',
        '--job-name=BDT_{0}'.format(i),
        '--time={0}'.format(time),
        '{0}/{1}/submitter_chunk_{2}.sh'.format(dt_string,folder,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







