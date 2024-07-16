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
inp =  "HOOK_PATH"  #f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/HOOK_DATE_TIME/HOOK_FILE_IN"
out =  f"/scratch/pahwagne/HOOK_DATE_TIME/HOOK_FILE_OUT"

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
nf = NanoFrame(inp, ) 

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
emptyFlag = nf["mu_pt"].any()

# sort the candidates according to:
# 1. we want opposite pi mu charge (priority)
# 2. we want opposite kk charge (2nd priority)
# 3. we want the dsMu mass to be physical (lower than Bs mass)
# 4. we sort after dsMu pt
  
# The logic goes as follows: if we have two candidates, one has the correct dsMu charge
# assignment but low dsMu pt, we still keep this event over wrong charge candidates with
# high pt, according to the priority list above. We do not throw away events! Worst case
# scenario is that we have only wrong charge candidates and take the one with highest pt
  
mySorting = lambda x : ( (x[iMuCharge]*x[iPiCharge] < 0.,x[iK1Charge]*x[iK2Charge] < 0. , x[iDsMuMass] < bsMass_, x[iDsMuPt]) )
 
  
#adapt these shapes to the nr of candidates in every event
nf["run"]             = np.array([np.array([run]  *len(nf[names[-1]][ev])) for ev,run  in enumerate(nf["run"]) ],             dtype=object)
nf["luminosityBlock"] = np.array([np.array([lumi] *len(nf[names[-1]][ev])) for ev,lumi in enumerate(nf["luminosityBlock"]) ], dtype=object)
nf["event"]           = np.array([np.array([evt]  *len(nf[names[-1]][ev])) for ev,evt  in enumerate(nf["event"]) ],           dtype=object)
nf["n"]               = np.array([np.array([n]    *len(nf[names[-1]][ev])) for ev,n    in enumerate(nf["n"]) ],               dtype=object)
try: nf["ngen"]       = np.array([np.array([n]    *len(nf[names[-1]][ev])) for ev,n    in enumerate(nf["ngen"]) ],            dtype=object)
except: pass
 

for i,name in enumerate(names):
  
  #Define all your branches
  nf[name] = nf[name][emptyFlag]
 
  # keep important indices for sorting
  if (name == "mu_charge"): iMuCharge = i    
  if (name == "pi_charge"): iPiCharge = i    
  if (name == "k1_charge"): iK1Charge = i    
  if (name == "k2_charge"): iK2Charge = i    
  if (name == "dsMu_m"   ): iDsMuMass = i    
  if (name == "dsMu_pt"  )  : iDsMuPt   = i    

nEntries = len(nf[names[-1]])
print(f"flattening {nEntries} events!")

# Lets rewrite the nfDir in the following shape and sort at the same time to avoid an extra loop:
#                   [        [ [this is one cand with all its variables pt, mass, m2miss, .... ]        #and we dump all candidates in an array]  sort this!                                # for each event]          
bestCands = np.array([ sorted([ [nf[name][ev][nCand] for name in names                       ] for nCand in range(len(nf[names[-1]][ev])) ], key = mySorting, reverse = True)[0] for ev in range(nEntries)])

# Lets sort them for each event (each element in allCands, as this is the outermost index (see above) and pick the best (index 0)
#bestCands = np.array([ listOfCands[0] for listOfCands in allCands ]) 
## For saving, we have to slice after the name rather than event
for i,name in enumerate(names):
  toSave[name] = bestCands[:,i]
  typesToSave[name] = type(toSave[name][0]) #take one element of to Save and get type
outfile = uproot.recreate(out)
outfile["tree"] = uproot.newtree(typesToSave)
outfile["tree"].extend(toSave)
 

