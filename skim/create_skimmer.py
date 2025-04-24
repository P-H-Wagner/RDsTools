import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import ROOT
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *



#############################################
# In this file we skim the files coming     #
# from the nanoAOD production. Additionally #
# we add the weighted reco branch.          #
#############################################

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def filesFromFolder(direc):
   filenames = os.listdir(direc)
   return [ direc + filename for filename in filenames ]

parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('channel')     # sig or hb
parser.add_argument('selection')   # pick selection
parser.add_argument('constrained', type = boolean_string)
parser.add_argument('prod')
args = parser.parse_args()

#python create_skimmer.py 24_04_2025_14_03_22 data base_wout_tv_25 True 25

queue = 'short' 
time = 10
nevents = -1
filesPerJob = 50

#naming and files

if args.constrained: 

  if args.prod == "24":
    #2024 production
    sig_cons = sig_cons_24
    if    args.channel == 'sig'    : naming = 'all_signals' 
    elif  args.channel == 'hb'     : naming = 'hb_inclusive'
    elif  args.channel == 'b0'     : naming = 'b0'          
    elif  args.channel == 'bs'     : naming = 'bs'          
    #elif  args.channel == 'lambdab': naming = 'lambdab'    
    elif  args.channel == 'bplus'  : naming = 'bplus'       
    else:                            naming = 'data'        

  else:
    #2025 production
    sig_cons = sig_cons_25
    if    args.channel == 'sig'    : naming = 'all_signals' 
    elif  args.channel == 'hb'     : naming = 'hb_inclusive'
    elif  args.channel == 'b0'     : naming = 'b0'          
    elif  args.channel == 'bs'     : naming = 'bs'          
    #elif  args.channel == 'lambdab': naming = 'lambdab'    
    elif  args.channel == 'bplus'  : naming = 'bplus'       
    else:                            naming = 'data'        

if not args.constrained:

  if    args.channel == 'sig'    : naming = 'all_signals' 
  elif  args.channel == 'hb'     : naming = 'hb_inclusive'
  elif  args.channel == 'b0'     : naming = 'b0'          
  elif  args.channel == 'bs'     : naming = 'bs'          
  #elif  args.channel == 'lambdab': naming = 'lambdab'    
  elif  args.channel == 'bplus'  : naming = 'bplus'       
  else:                            naming = 'data'        


if args.channel == 'fakes':
  naming = 'fakes'; files = fakes;

# get chain to access variable names (need all files bc sometimes they are empty for data)
directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/*"
chain = ROOT.TChain("tree")
chain.Add(directory)

# get the inputfiles
inputfiles = filesFromFolder(directory[:-1]) # remove * used for the chain

#define weights for math solution (take the most actual signal sample)
signals = ROOT.TChain("tree")

if args.constrained:       input_sig = sig_cons
if not args.constrained:   input_sig = sig_unc
if args.channel == 'fakes': input_sig = sig_cons #fakes so far only produced for constrained

for sig in input_sig:
  signal_directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{sig}/*"
  signals.Add(signal_directory)

# get all non-nan solutions
def getWeights(name):
  print("checking", name)
  ntot   = signals.GetEntries("(gen_sig == 0)")                                                                          #check that reco is not nan!
  nReco1 = signals.GetEntries(f"(abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && (gen_sig == 0) && !({name}_reco_1 != {name}_reco_1)")
  w1 = nReco1/ntot
  print("w1 is: ", w1)
  return w1 , 1 - w1


#branches and names
branches = chain.GetListOfBranches()
names    = [branch.GetName() for branch in branches]

# errors and logs

if not os.path.exists(f"./{args.date_time}_{args.selection}/"):
  os.makedirs(f"{args.date_time}_{args.selection}/logs")
  os.makedirs(f"{args.date_time}_{args.selection}/errs")

#write for loop as string to print into the skimmer file

forLoop = "df \\"


w1, w2 = getWeights("bs_pt")

for name in names:
  #check if reco is available
  if "reco_1" in name:
 
    # reco_2 string:
    reco_2 = name.replace("_reco_1", "_reco_2")
    reco_w = name.replace("_reco_1", "_reco_weighted")
 
    #core = name[:-7] # get var name without "_reco_1"
    #print("core name is: ", core, "var name is: ", name)
    #w1, w2 = getWeights(core)
    #forLoop += f"\n.Define(\"{core}_reco_weighted\",\"{w1}*{core}_reco_1 + {w2}*{core}_reco_2\")\\"  
    forLoop += f"\n.Define(\"{reco_w}\",\"{w1}*{name} + {w2}*{reco_2}\")\\"  

forLoop += f"\n.Filter(selec).Snapshot(\"tree\", destination)"
print(forLoop)

# loop over all to be skimmed files

for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_skimmer.py", "rt")
  #file to write to
  cfg = open(f"{args.date_time}_{args.selection}/skimmer_{i}.py","wt")

  fout = f"/scratch/pahwagne/skimmed_{args.date_time}_{args.selection}/skimmed_{args.date_time}_{args.selection}_chunk_{i}.root"
  

  for line in temp:
    if   "HOOK_FILE_IN"  in line: 
      for f in fin:
        cfg.write(line.replace("HOOK_FILE_IN", f))
      continue

    elif "HOOK_DATE_TIME"   in line: line = line.replace("HOOK_DATE_TIME", args.date_time)
    elif "HOOK_SELECTION"   in line: line = line.replace("HOOK_SELECTION", baselines[args.selection])
    elif "HOOK_NEW_BRANCH"  in line: line = line.replace('"HOOK_NEW_BRANCH"', forLoop)
    elif "HOOK_FILE_OUT"    in line: line = line.replace("HOOK_FILE_OUT", fout)
  
    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'eval "$(conda shell.bash hook)"',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/skimmed_{args.date_time}_{args.selection}',
         f'python /work/pahwagne/RDsTools/skim/{args.date_time}_{args.selection}/skimmer_{i}.py',
         f'xrdcp {fout} root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.date_time}/skimmed_{args.selection}_{args.date_time}_chunk_{i}.root',
         f'rm -r {fout}',
         '',
     ])
  
  with open(f"{args.date_time}_{args.selection}/skimmer_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)
  
  
  command_sh_batch = ' '.join([
  
        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o {args.date_time}_{args.selection}/logs/log_{i}.txt',
        f'-e {args.date_time}_{args.selection}/errs/err_{i}.txt',
        #'--mem=1500M',
        f'--job-name=SKIM_{args.channel}',
        #f'--time={time}',
        f'{args.date_time}_{args.selection}/skimmer_{i}.sh',
     ])
  
  print(command_sh_batch)
  os.system(command_sh_batch)







