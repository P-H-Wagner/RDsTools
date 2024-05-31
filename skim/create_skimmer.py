import os
from glob import glob
from datetime import datetime
import sys 
import argparse
import re #for pattern matchin
import ROOT


parser = argparse.ArgumentParser()

parser.add_argument('date_time')
parser.add_argument('channel')   # sig or hb
parser.add_argument('selection') # pick selection
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

# get cahin to access variable names (need all files bc sometimes they are empty for data)
directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.date_time}/*"
chain = ROOT.TChain("tree")
chain.Add(directory)

#define weights for math solution
signal_directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/13_05_2024_14_39_06/*"
signals = ROOT.TChain("tree")
signals.Add(signal_directory)

# get all non-nan solutions
def getWeights(name):
  print("checking", name)
  ntot   = signals.GetEntries("(gen_sig == 0)")
  nReco1 = signals.GetEntries(f"(abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && (gen_sig == 0)")
  w1 = nReco1/ntot
  print("w1 is: ", w1)
  return w1 , 1 - w1


#branches and names
branches = chain.GetListOfBranches()
names    = [branch.GetName() for branch in branches]

# errors and logs
os.makedirs(f"{args.date_time}/logs")
os.makedirs(f"{args.date_time}/errs")

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
cfg = open(f"{args.date_time}/skimmer.py","wt")

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
       f'mkdir -p /scratch/pahwagne',
       f'python /work/pahwagne/RDsTools/comb/{args.date_time}/skimmer.py',
       f'xrdcp /scratch/pahwagne/skimmed_{args.selection}_{args.date_time}.root root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.date_time}/skimmed_{args.selection}_{args.date_time}.root',
       f'rm -r /scratch/pahwagne',
       '',
   ])

with open(f"{args.date_time}/skimmer.sh", "wt") as flauncher:
  flauncher.write(to_write)


command_sh_batch = ' '.join([

      'sbatch',
      f'-p {queue}',
      '--account=t3',
      f'-o {args.date_time}/logs/log.txt',
      f'-e {args.date_time}/errs/err.txt',
      #'--mem=1500M',
      f'--job-name=SKIM_{args.channel}',
      #f'--time={time}',
      f'{args.date_time}/skimmer.sh',
   ])

print(command_sh_batch)
os.system(command_sh_batch)







