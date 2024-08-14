import os
from glob import glob
from datetime import datetime
import sys
import argparse
import re #for pattern matchin
import ROOT
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

parser = argparse.ArgumentParser()
parser.add_argument('selection')
args = parser.parse_args()
sel = args.selection


#skim all!

def run(cmd):
  print(cmd)
  os.system(cmd)

##unconstrained
for isig in sig_unc:
  cmd = f"python create_skimmer.py {isig} sig {sel} False"
  run(cmd)
for isig in hb_unc:
  cmd = f"python create_skimmer.py {isig} hb {sel} False"
  run(cmd)
for isig in bs_unc:
  cmd = f"python create_skimmer.py {isig} bs {sel} False"
  run(cmd)
for isig in b0_unc:
  cmd = f"python create_skimmer.py {isig} b0 {sel} False"
  run(cmd)
for isig in bplus_unc:
  cmd = f"python create_skimmer.py {isig} bplus {sel} False"
  run(cmd)
for isig in data_unc:
  cmd = f"python create_skimmer.py {isig} data {sel} False"
  run(cmd)

for isig in sig_cons:
  cmd = f"python create_skimmer.py {isig} sig {sel} True"
  run(cmd)
for isig in hb_cons:
  cmd = f"python create_skimmer.py {isig} hb {sel} True"
  run(cmd)
for isig in bs_cons:
  cmd = f"python create_skimmer.py {isig} bs {sel} True"
  run(cmd)
for isig in b0_cons:
  cmd = f"python create_skimmer.py {isig} b0 {sel} True"
  run(cmd)
for isig in bplus_cons:
  cmd = f"python create_skimmer.py {isig} bplus {sel} True"
  run(cmd)
for isig in data_cons:
  cmd = f"python create_skimmer.py {isig} data {sel} True"
  run(cmd)

