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

parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('channel')     # sig or hb
parser.add_argument('selection')   # pick selection
parser.add_argument('constrained', type = boolean_string)
args = parser.parse_args()

queue = 'short' 
time = 10
nevents = -1

#naming and files

if args.constrained: 

  #file_data  = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in data_cons ]
  #file_sig   = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in sig_cons  ]
  #file_hb    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in hb_cons   ]
  #file_b0    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in b0_cons   ]
  #file_bs    = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bs_cons   ]
  #file_bplus = [f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{f}/skimmed_base_wout_tv_{f}.root" for f in bplus_cons]


  if    args.channel == 'sig'    : naming = 'all_signals' ; files = sig_cons;
  elif  args.channel == 'hb'     : naming = 'hb_inclusive'; files = hb_cons;
  elif  args.channel == 'b0'     : naming = 'b0'          ; files = b0_cons;
  elif  args.channel == 'bs'     : naming = 'bs'          ; files = bs_cons;
  #elif  args.channel == 'lambdab': naming = 'lambdab'     ; files = lamdbab_cons;
  elif  args.channel == 'bplus'  : naming = 'bplus'       ; files = bplus_cons;
  else:                            naming = 'data'        ; files = data_cons;

if not args.constrained:

  if    args.channel == 'sig'    : naming = 'all_signals' ; files = sig_unc;
  elif  args.channel == 'hb'     : naming = 'hb_inclusive'; files = hb_unc;
  elif  args.channel == 'b0'     : naming = 'b0'          ; files = b0_unc;
  elif  args.channel == 'bs'     : naming = 'bs'          ; files = bs_unc;
  #elif  args.channel == 'lambdab': naming = 'lambdab'     ; files = lamdbab_unc;
  elif  args.channel == 'bplus'  : naming = 'bplus'       ; files = bplus_unc;
  else:                            naming = 'data'        ; files = data_unc;



if args.channel == 'fakes':
  naming = 'fakes'; files = fakes;


# get cahin to access variable names (need all files bc sometimes they are empty for data)
directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/*"
chain = ROOT.TChain("tree")
chain.Add(directory)

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
os.makedirs(f"{args.date_time}_{args.selection}/logs")
os.makedirs(f"{args.date_time}_{args.selection}/errs")

#write for loop as string to print into the skimmer file

forLoop = "df \\"


w1, w2 = getWeights("bs_pt")

for name in names:

  #check if reco is available
  if "reco_1" in name:
    core = name[:-7] # get var name without "_reco_1"
    #w1, w2 = getWeights(core)
    forLoop += f"\n.Define(\"{core}_reco_weighted\",\"{w1}*{core}_reco_1 + {w2}*{core}_reco_2\")\\"  

forLoop += f"\n.Filter({args.selection}).Snapshot(\"tree\", destination)"
print(forLoop)

#template
temp = open("temp_skimmer.py", "rt")
#file to write to
cfg = open(f"{args.date_time}_{args.selection}/skimmer.py","wt")

for line in temp:
  if "HOOK_DATE_TIME"   in line: line = line.replace("HOOK_DATE_TIME", args.date_time)
  if "HOOK_SELECTION"   in line: line = line.replace("HOOK_SELECTION", args.selection)
  if "HOOK_NEW_BRANCH"  in line: line = line.replace('"HOOK_NEW_BRANCH"', forLoop)

  cfg.write(line)

temp.close()
cfg.close()

to_write = '\n'.join([
       '#!/bin/bash',
       'eval "$(conda shell.bash hook)"',
       'conda activate /work/pahwagne/environments/hammer3p8',
       f'mkdir -p /scratch/pahwagne/skimmed_{args.date_time}_{args.selection}',
       f'python /work/pahwagne/RDsTools/skim/{args.date_time}_{args.selection}/skimmer.py',
       f'xrdcp /scratch/pahwagne/skimmed_{args.selection}_{args.date_time}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.date_time}/skimmed_{args.selection}_{args.date_time}.root',
       f'rm -r /scratch/pahwagne/skimmed_{args.date_time}_{args.selection}',
       '',
   ])

with open(f"{args.date_time}_{args.selection}/skimmer.sh", "wt") as flauncher:
  flauncher.write(to_write)


command_sh_batch = ' '.join([

      'sbatch',
      f'-p {queue}',
      '--account=t3',
      f'-o {args.date_time}_{args.selection}/logs/log.txt',
      f'-e {args.date_time}_{args.selection}/errs/err.txt',
      #'--mem=1500M',
      f'--job-name=SKIM_{args.channel}',
      #f'--time={time}',
      f'{args.date_time}_{args.selection}/skimmer.sh',
   ])

print(command_sh_batch)
os.system(command_sh_batch)







