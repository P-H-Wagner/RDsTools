import os
import numpy as np
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
# we add the weighted reco branch and the   #
# pileup weights.                           #
#############################################

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

def filesFromFolder(direc):
   filenames = os.listdir(direc)
   return [ direc + filename for filename in filenames ]

parser = argparse.ArgumentParser()

parser.add_argument("--date_time", required = True)
parser.add_argument("--channel"  , required = True)  
parser.add_argument("--selection", required = True)  
parser.add_argument("--hammered" , required = True, type = boolean_string)   

args = parser.parse_args()

#python create_skimmer.py --date_time=20260227_105644 --channel=data --selection=base_wout_tv_25 --hammered=False
#python create_skimmer.py --date_time=20_03_2026_09_37_22 --channel=sig --selection=minimal --hammered=True

####################

# datetime logic 
dt_in = args.date_time

now    = datetime.now()
dt     = now.strftime("%d_%m_%Y_%H_%M_%S") 

####################

if args.hammered:
  dt_out = dt #save as a new dir without overwriting the old hammer
else:
  dt_out = dt_in 

queue = 'short' 
time = 20
nevents = -1
filesPerJob = 10

#naming and files
if    args.channel == 'sig'    : naming = 'all_signals' 
elif  args.channel == 'hb'     : naming = 'hb_inclusive'
elif  args.channel == 'b0'     : naming = 'b0'          
elif  args.channel == 'bs'     : naming = 'bs'          
#elif  args.channel == 'lambdab': naming = 'lambdab'    
elif  args.channel == 'bplus'  : naming = 'bplus'       
else:                            naming = 'data'        

# get chain to access variable names (need all files bc sometimes they are empty for data)
if args.hammered:
  directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/signal_default_{dt_in}/*"
else:
  directory = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dt_in}/*"

chain = ROOT.TChain("tree")
chain.Add(directory)

# get the inputfiles
inputfiles = filesFromFolder(directory[:-1]) # remove * used for the chain
#inputfiles = inputfiles[:1]

#branches and names
branches = chain.GetListOfBranches()
names    = [branch.GetName() for branch in branches]

# errors and logs

if not os.path.exists(f"./{dt_in}_{args.selection}/"):
  os.makedirs(f"{dt_in}_{args.selection}/logs")
  os.makedirs(f"{dt_in}_{args.selection}/errs")


###############################
# Corrections and systematics #
###############################

if args.channel == "sig":

  #rename time branch to avoid clashes with cpp time library ...
   
  forLoop = f"df = df.Define  ('bs_time'    ,'get_time(gen_pv_x, gen_pv_y, gen_pv_z, gen_sv_x, gen_sv_y, gen_sv_z, gen_bs_boost)') \n"
  forLoop += f"df \\"

  #todo add up/down for PU
  forLoop += f"\n.Define  ('w_pu'              ,'get_w_pu              (gen_numPUInts)           ') \\"
  forLoop += f"\n.Define  ('w_pu_up'           ,'get_w_pu_up           (gen_numPUInts)           ') \\"
  forLoop += f"\n.Define  ('w_pu_down'         ,'get_w_pu_down         (gen_numPUInts)           ') \\"

  #forLoop += f"\n.Redefine('trigger_sf'       ,'get_trigger_sf    (mu_pt, dxy_mu_sig_pv)   ') \\"
  forLoop += f"\n.Define  ('trigger_sf'        ,'get_trigger_sf    (mu_pt, dxy_mu_sig_pv)   ') \\"
  forLoop += f"\n.Define  ('trigger_sf_err_hi' ,'get_trigger_sf_err_hi(mu_pt, dxy_mu_sig_pv)   ') \\"
  forLoop += f"\n.Define  ('trigger_sf_err_lo' ,'get_trigger_sf_err_lo(mu_pt, dxy_mu_sig_pv)   ') \\"
  #forLoop += f"\n.Define  ('trigger_sf_err'   ,'get_trigger_sf_err(mu_pt, dxy_mu_sig_pv)   ') \\"
  forLoop += f"\n.Define  ('w_bs_tau'          ,'get_w_bs_tau      (bs_time)                ') \\"
  forLoop += f"\n.Define  ('w_bs_tau_up'       ,'get_w_bs_tau_up   (bs_time)                ') \\"
  forLoop += f"\n.Define  ('w_bs_tau_down'     ,'get_w_bs_tau_down (bs_time)                ') \\"

if args.channel == "hb":
  forLoop  = f"df.Define  ('trigger_sf'        ,'get_trigger_sf    (mu_pt, dxy_mu_sig_pv)') \\"
  forLoop += f"\n.Define  ('w_pu'              ,'get_w_pu              (gen_numPUInts)           ') \\"
  forLoop += f"\n.Define  ('w_pu_up'           ,'get_w_pu_up           (gen_numPUInts)           ') \\"
  forLoop += f"\n.Define  ('w_pu_down'         ,'get_w_pu_down         (gen_numPUInts)           ') \\"
  forLoop += f"\n.Define  ('trigger_sf_err_hi' ,'get_trigger_sf_err_hi(mu_pt, dxy_mu_sig_pv)') \\"
  forLoop += f"\n.Define  ('trigger_sf_err_lo' ,'get_trigger_sf_err_lo(mu_pt, dxy_mu_sig_pv)') \\"
  #forLoop += f"\n.Define  ('trigger_sf_err'   ,'get_trigger_sf_err(mu_pt, dxy_mu_sig_pv)') \\"

if args.channel == "data":
  forLoop = "df \\"
 
  if (dt_in == "20260227_105644") or (dt_in == "20260227_104615") or (dt_in == "20260227_102607"):
    #dirty fix
    #since the period D was produced on crb earlier this year without the 8p5 and 10p5 parts, lets add a column manually
    #we set the trigger zero all the time, since it was never on in period D!
    
    forLoop += f"\n.Define  ('mu8p5_3p5'  , '0.0') \\"
    forLoop += f"\n.Define  ('mu10p5_3p5' , '0.0') \\"

###############################
# take snapshot               #
###############################

forLoop += f"\n.Filter(selec).Snapshot(\"tree\", destination)"
print(forLoop)

# loop over all to be skimmed files

for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
 
  #template
  temp = open("temp_skimmer.py", "rt")
  #file to write to
  cfg = open(f"{dt_in}_{args.selection}/skimmer_{i}.py","wt")


  fout    = f"/scratch/pahwagne/skimmed_{dt_out}_{args.selection}/skimmed_{dt_out}_{args.selection}_chunk_{i}.root"
 
  if args.hammered:
    pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/signal_default_{dt_out}/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/signal_default_{dt_out}_D1_sf_only/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/signal_default_{dt_out}_A1_sf_only/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/signal_default_{dt_out}_D1_sf_only/")
    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/hammer/{args.selection}/signal_default_{dt_out}/")
  else: 
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_D1_sf_only/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_ricc_D_only/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_fullBPark/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"
    #pnfs_loc = f"root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_A1_sf_only/skimmed_{args.selection}_{dt_out}_chunk_{i}.root"

    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}/")
    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_D1_sf_only/")
    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_ricc_D_only/")
    os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_fullBPark/")
    #os.system(f"mkdir -p /pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{args.selection}/{dt_out}_A1_sf_only/")


  for line in temp:
    if   "HOOK_FILE_IN"  in line: 
      for f in fin:
        cfg.write(line.replace("HOOK_FILE_IN", f))
      continue

    elif "HOOK_DATE_TIME"   in line: line = line.replace("HOOK_DATE_TIME"   , dt_new)
    elif "HOOK_SELECTION"   in line: line = line.replace("HOOK_SELECTION"   , baselines[args.selection])
    elif "HOOK_NEW_BRANCH"  in line: line = line.replace('"HOOK_NEW_BRANCH"', forLoop)
    elif "HOOK_FILE_OUT"    in line: line = line.replace("HOOK_FILE_OUT"    , fout)
  
    cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         'eval "$(conda shell.bash hook)"',
         #'conda activate /work/pahwagne/environments/skim_env',
         'conda activate /work/pahwagne/environments/hammer3p8',
         f'mkdir -p /scratch/pahwagne/skimmed_{dt_out}_{args.selection}',
         f'python /work/pahwagne/RDsTools/skim/{dt_in}_{args.selection}/skimmer_{i}.py',
         f'xrdcp {fout} {pnfs_loc}',
         f'rm -r {fout}',
         '',
     ])
  
  with open(f"{dt_in}_{args.selection}/skimmer_{i}.sh", "wt") as flauncher:
    flauncher.write(to_write)
  
  
  command_sh_batch = ' '.join([
  
        'sbatch',
        f'-p {queue}',
        '--account=t3',
        f'-o {dt_in}_{args.selection}/logs/log_{i}.txt',
        f'-e {dt_in}_{args.selection}/errs/err_{i}.txt',
        #'--mem=1500M',
        f'--job-name=SKIM_{args.channel}',
        f'--time={time}',
        f'{dt_in}_{args.selection}/skimmer_{i}.sh',
     ])
  
  print(command_sh_batch)
  os.system(command_sh_batch)







