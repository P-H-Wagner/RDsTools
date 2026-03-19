import ROOT
import argparse
import os
import sys

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import * 

import numpy as np

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('constrained', type=boolean_string) #use constrained sampels?
args = parser.parse_args()

if args.constrained: data = data_cons
else: data = data_unc

sample_list = []
for f in data:
  samples = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{f}/*"
  sample_list += glob.glob(samples)

#branches to load
branches = [
"run",
"luminosityBlock",
]


for i,sample in enumerate(sample_list):

  with uproot.open(sample) as f:
    
      tree = f["tree"]

      if i == 0:
        #starting value
        df   = tree.arrays(branches, library = "pd", entry_start=None, entry_stop=None, cut=selection[chan]) 

      else:
        #append
        df = df.append(tree.arrays(branches, library = "pd", entry_start=None, entry_stop=None, cut=selection[chan]) ) 


