import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import subprocess

#example task: 240802_091822:pahwagne_crab_20240802_111807
#example command line: python create_submitter.py -f 0 -i  240802_091822 -p 2 20240802_111807 data t2
#example command line: python create_submitter.py -f 0 -i  250217_105711 -p 1 20250217_114156 data t2
#example command line: python create_submitter.py 16_03_2025_20_56_34 sig t3
parser = argparse.ArgumentParser()

parser.add_argument('date_time') # second number sequence of crab task name
parser.add_argument('channel') # sig, data, hb, bs , ...
parser.add_argument('site')  # t3 or t2 (from crab)
parser.add_argument('-f','--folder')  # nr of directory on the t2, f.e. nDir = 2 means: 0002 
parser.add_argument('-i','--crab_id') # first number sequence of crab task name 
parser.add_argument('-p','--part')    # 1,2,3,4,5 - specifies the part of the dataset (BPH1, BPH2, ..., BPH5) OR the binning for qcd (e.g. 50to100)
args = parser.parse_args()


## experience: if crab file aleady contains 10 files per job, this is the max the flattening can handle without setting the memory up to >2GB


queue = 'short' 
#time = 60
nevents = -1

#naming
if args.channel == 'sig': naming = 'all_signals'
elif args.channel == 'hb': naming = 'hb_inclusive'
elif args.channel == 'b0': naming = 'b0'
elif args.channel == 'bs': naming = 'bs'
elif args.channel == 'lambdab': naming = 'lambdab'
elif args.channel == 'fakes': naming = 'fakes'
elif args.channel == 'bplus': naming = 'bplus'
elif args.channel == 'dsmu': naming = 'dsmu'
elif args.channel == 'dsmu_isgw2': naming = 'dsmu_isgw2'
elif args.channel == 'dstau': naming = 'dstau'
elif args.channel == 'dsstarmu': naming = 'dsstarmu'
elif args.channel == 'dsstarmu_isgw2': naming = 'dsstarmu_isgw2'
elif args.channel == 'dsstartau': naming = 'dsstartau'
else: naming = 'data'


template   = "temp_flatNanoProducer.py"

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

#load files
gen = False

if ((args.channel == "dstau") or (args.channel == "dsstartau") or (args.channel == "dsmu") or (args.channel == "dsstarmu") or (args.channel == "dsmu_isgw2") or (args.channel == "hb_gen") or (args.channel == "dsstarmu_isgw2")):
  path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/inspector/" 
  directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/inspector/{args.date_time}/"
  inputfiles = os.listdir(directory)
  template   = "temp_flatNanoProducer_gen.py"
  gen = True
  log_path = f"{args.date_time}"   

elif args.site == 't3' and not gen:
  path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/" 
  directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/{args.date_time}/"
  inputfiles = os.listdir(directory)

  log_path = f"{args.date_time}"   

else:

  # Convert Date-Time into crab naming:
  
  months = {"02": "Feb", "06": "Jun", "07": "Jul", "08": "Aug", "09": "Sep", "10": "Oct"} # add more if needed
  
  year  = args.date_time[0:4]
  month = months[args.date_time[4:6]]
  day   = args.date_time[6:8]
  
  crab_date = year+month+day

  try: folder = args.folder
  except: print("Please give the folder number on the T2 storage for this task!")

  try: crab_id = args.crab_id
  except: print("Please give also the crab id next to date time for this task!")

  try: part = args.part
  except: print("Please give also the data part for this task!")

  if args.channel == "fakes":
    path        = f"root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat//store/user/pahwagne/{crab_date}_{part}/QCD_HT{part}_TuneCP5_13TeV-madgraphMLM-pythia8/crab_{args.date_time}/{args.crab_id}/000{folder}/"
  else:
    path        = f"root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat//store/user/pahwagne/{crab_date}/ParkingBPH{part}/crab_{args.date_time}/{args.crab_id}/000{folder}/"

  files       = list_files_gfal(path)
  inputfiles  = [path + f for f in files] #inlcude the path in the list :-)
 
  #path       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/"
  #directory  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/crab/{args.date_time}/"
  #inputfiles = os.listdir(directory)

  log_path = f"{args.date_time}_{args.folder}"   


#inputfiles = inputfiles[0:1] #debug
#print(inputfiles)

if not os.path.exists(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}"):

  os.makedirs(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}")

if not os.path.exists(f"flat_{log_path}"):
  os.makedirs(f"flat_{log_path}/logs")
  os.makedirs(f"flat_{log_path}/errs")

#inputfiles = inputfiles[0:5]
#print(inputfiles)
for fin in inputfiles:
  print("here")

  if args.site == 't3':
    #get chunk nr
    pattern = r'(.+)_chunk_(\d+)\.root' 
    match = re.search(pattern, fin)
    if match: chunkNr = match.group(2)

  else:
    #get chunk nr
    pattern = r'([^_]+)_crab_(\d+)\.root'
    match = re.search(pattern, fin)
    if match: chunkNr = match.group(2)      

  print(chunkNr)

  #template
  temp = open(template, "rt")
  #file to write to
  cfg = open(f"flat_{log_path}/cfg_chunk_{chunkNr}.py","wt")
  #file to save things
  fout = f"{naming}_flatChunk_{chunkNr}.root"  

  for line in temp:
    if "HOOK_DATE_TIME" in line: line = line.replace("HOOK_DATE_TIME", args.date_time)
    if "HOOK_LOG_PATH" in line: line  = line.replace("HOOK_LOG_PATH", log_path)
    if "HOOK_FILE_IN" in line:   line = line.replace("HOOK_FILE_IN", fin)
    if "HOOK_FILE_OUT" in line:  line = line.replace("HOOK_FILE_OUT", fout)
    if "HOOK_CHUNK" in line:     line = line.replace("HOOK_CHUNK", chunkNr)
    if "HOOK_CHANNEL" in line:   line = line.replace("HOOK_CHANNEL", naming)
    if "HOOK_PATH" in line:      
      if args.site == "t3": line = line.replace("HOOK_PATH", path + '/' + args.date_time + '/' + fin)
      else                : line = line.replace("HOOK_PATH", fin)

    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'cd /work/pahwagne/RDsTools/flatNano',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'cp nanoframe.py flat_{log_path}/nanoframe.py',
         f'mkdir -p /scratch/pahwagne/{log_path}',
         f'python flat_{log_path}/cfg_chunk_{chunkNr}.py',
         f'xrdcp /scratch/pahwagne/{log_path}/{fout} root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/{fout}',
         f'rm /scratch/pahwagne/{log_path}/{fout}',
         '',
     ])

  with open(f"flat_{log_path}/submitter_chunk_{chunkNr}.sh", "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o flat_{log_path}/logs/chunk_{chunkNr}.log',
        f'-e flat_{log_path}/errs/chunk_{chunkNr}.err',
        f'--mem=3000M',
        f'--job-name=FLAT_{chunkNr}_{naming}',
        #f'--time=60',
        f'flat_{log_path}/submitter_chunk_{chunkNr}.sh',
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







