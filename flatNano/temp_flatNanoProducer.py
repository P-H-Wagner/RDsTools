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

#Single File
#fileName = "/work/pahwagne/test/multi_instances.root" 
#load it into root file to get all branch names
#rootFile = ROOT.TFile.Open(fileName)
#rootTree = rootFile.Get("Events")

#before we put the bs reconstruction methods into functions
#inp = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/11_02_2024_21_56_13/"
#and after


# the input can be a single file (t3 production) or a list of files (crab)
#inp =  "HOOK_PATH"  #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/HOOK_DATE_TIME/HOOK_FILE_IN"
inp =  f"/scratch/pahwagne/HOOK_LOG_PATH/HOOK_FILE_SRC"
out =  f"/scratch/pahwagne/HOOK_LOG_PATH/HOOK_FILE_OUT"
#load it into root file to get all branch names
rootFile = ROOT.TFile.Open(inp)
rootTree = rootFile.Get("Events")


# nr of events
#nEntries = rootTree.GetEntries()
#print(f" ==> Start flatteing {nEntries} Events of chunk HOOK_CHUNK")
  
##################
### Flattening  ##
##################
  
#will hold the nanoAOD values -> can handle jagged arrays!!
nf_bare = NanoFrame(inp, ) 

#dictionary which holds everything
#nfDir = {}
#dictionary which holds final candidate and their types
toSave = {}
typesToSave = {}

  
#branches and names
branches = rootTree.GetListOfBranches()
names    = [branch.GetName() for branch in branches]
# remove the last few branches (f.e. fixedGridRhoFastjetAll) 
# Do we need them ? otherwise I have to adapt the shape as below
names = [name for name in names if not (name.startswith("fixed") or name == "nmuon" or name.startswith("muon")) ] 


#Bs mass
bsMass_ = 5.36688

#chekc if not empty
emptyFlag = nf_bare["mu_pt"].any()

# sort the candidates according to:
# 1. we want opposite pi mu charge (priority)
# 2. we want opposite kk charge (2nd priority)
# 3. we want the dsMu mass to be physical (lower than Bs mass)
# 4. we sort after dsMu pt
  
# The logic goes as follows: if we have two candidates, one has the correct dsMu charge
# assignment but low dsMu pt, we still keep this event over wrong charge candidates with
# high pt, according to the priority list above. We do not throw away events! Worst case
# scenario is that we have only wrong charge candidates and take the one with highest pt
  
#mySorting = lambda x : ( (x[iMuCharge]*x[iPiCharge] < 0.,x[iK1Charge]*x[iK2Charge] < 0. , x[iDsMuMass] < bsMass_, x[iDsMuPt]) )
#mySorting = lambda x : ( (x[iMuCharge]*x[iPiCharge] < 0.,x[iK1Charge]*x[iK2Charge] < 0. , x[iDsMuPt]) )
mySorting = lambda x : ( ( x[iDsMuPt]) )
 
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


reference = nf[names[-1]]

for key in ["run", "luminosityBlock", "event", "n"]:
    #nf[key] = ak.broadcast_arrays(nf[key], reference)[0]
    nf[key] = np.array([np.array([val]  *len(nf[names[-1]][ev])) for ev,val  in enumerate(nf[key]) ],             dtype=object)

if "ngen" in nf.keys:
    #nf["ngen"] = ak.broadcast_arrays(nf["ngen"], reference)[0]
    nf[key] = np.array([np.array([val]  *len(nf[names[-1]][ev])) for ev,val  in enumerate(nf[key]) ],             dtype=object)

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
 

