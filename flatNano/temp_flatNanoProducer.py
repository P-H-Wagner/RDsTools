import awkward as ak
import numpy as np
import uproot
from nanoframe import NanoFrame
import os
import particle
import pandas as pd
import uproot_methods
import ROOT
from pdb import set_trace
from root_pandas import to_root
from uproot_methods import TLorentzVectorArray
from uproot_methods import TVector3Array
from uproot_methods import TLorentzVector
from uproot_methods import TVector3
from scipy.constants import c as speed_of_light
import uproot
import math
from array import array
import os
import sys


# the input can be a single file (t3 production) or a list of files (crab)
#inp =  "HOOK_PATH"  #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/HOOK_DATE_TIME/HOOK_FILE_IN"
inp =  f"/scratch/pahwagne/HOOK_LOG_PATH/HOOK_FILE_SRC"
out =  f"/scratch/pahwagne/HOOK_LOG_PATH/HOOK_FILE_OUT"
#load it into root file to get all branch names
rootFile = ROOT.TFile.Open(inp)
rootTree = rootFile.Get("Events")

##################
### Flattening  ##
##################
  
#will hold the nanoAOD values -> can handle jagged arrays!!
nf_bare = NanoFrame(inp, ) 

#dictionary which holds final candidate and their types
toSave = {}
typesToSave = {}

#branches and names
branches = rootTree.GetListOfBranches()
names    = [branch.GetName() for branch in branches]
# remove some branches (f.e. fixedGridRhoFastjetAll) 
names = [name for name in names if not (name.startswith("fixed") or name == "nmuon" or name.startswith("muon")) ] 

#Bs mass
bsMass_ = 5.36688

#chekc if not empty
emptyFlag = nf_bare["mu_pt"].any()

###############
# SIGN FILTER #
###############

exclude = ["run","luminosityBlock","event","n","ngen"]


if "HOOK_SIGN" == "WS":
  mask = ((nf_bare.array("pi_charge") * nf_bare.array("mu_charge") > 0) & (nf_bare.array("k1_charge") * nf_bare.array("k2_charge") < 0))
  dummy = {key: nf_bare.array(key)[mask] for key in names if key not in exclude}

elif "HOOK_SIGN" == "NONE":
  dummy = {key: nf_bare.array(key)       for key in names if key not in exclude}

else:
  mask = ((nf_bare.array("pi_charge") * nf_bare.array("mu_charge") < 0) & (nf_bare.array("k1_charge") * nf_bare.array("k2_charge") < 0))
  dummy = {key: nf_bare.array(key)[mask] for key in names if key not in exclude}

ws_any_mask = ((nf_bare.array("pi_charge") * nf_bare.array("mu_charge") > 0) | (nf_bare.array("k1_charge") * nf_bare.array("k2_charge") > 0))


dummy["run"]             = nf_bare.array("run")
dummy["luminosityBlock"] = nf_bare.array("luminosityBlock")
dummy["event"]           = nf_bare.array("event")
dummy["n"]               = nf_bare.array("n")
try: dummy["ngen"]       = nf_bare.array("ngen")
except: pass
nf      = NanoFrame(dummy)

nf["run"]             = ak.JaggedArray.fromiter([[run] * len(nf[names[-1]][ev]) for ev, run in enumerate(nf["run"])])
nf["luminosityBlock"] = ak.JaggedArray.fromiter([[luminosityBlock] * len(nf[names[-1]][ev]) for ev, luminosityBlock in enumerate(nf["luminosityBlock"])])
nf["event"]           = ak.JaggedArray.fromiter([[event] * len(nf[names[-1]][ev]) for ev, event in enumerate(nf["event"])])
nf["n"]               = ak.JaggedArray.fromiter([[n] * len(nf[names[-1]][ev]) for ev, n in enumerate(nf["n"])])
try:
  nf["ngen"] = ak.JaggedArray.fromiter([[n] * len(nf[names[-1]][ev]) for ev, n in enumerate(nf["ngen"])])
  pass

emptyFlag = nf["mu_pt"].any() 

for i,name in enumerate(names):
  
  #Define all your branches
  nf[name] = nf[name][emptyFlag]
 
  # keep important indices for sorting
  if (name == "mu_charge"): iMuCharge   = i    
  if (name == "pi_charge"): iPiCharge   = i    
  if (name == "k1_charge"): iK1Charge   = i    
  if (name == "k2_charge"): iK2Charge   = i    
  if (name == "dsMu_m"   ): iDsMuMass   = i    
  if (name == "dsMu_pt"  ): iDsMuPt     = i    

nEntries = len(nf[names[-1]])
print(f"flattening {nEntries} events!")

## For saving, we have to slice after the name rather than event
for i,name in enumerate(names):

  # Flatten completely --> keep all candidates per event
  nf[name] = nf[name].flatten()
 
  # unsigned types raise errors 
  if nf[name].dtype.kind == "u":
      nf[name] = nf[name].astype(np.int64)
  
  # Save the dataframe as a dictionary
  toSave[name] = nf[name]
  # save the type
  typesToSave[name] = type(toSave[name][0])


outfile = uproot.recreate(out)
outfile["tree"] = uproot.newtree(typesToSave)
outfile["tree"].extend(toSave)
 

