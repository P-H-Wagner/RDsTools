#from nanoAOD root files, to flat ntuples
#from coffea.analysis_objects import JaggedCandidateArray
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

#test input with multple indices
fileName = "/work/pahwagne/test/multi_instances.root" 

#load it into root file to get all branch names
rootFile = ROOT.TFile.Open(fileName)
rootTree = rootFile.Get("Events")

# nr of events
nEntries = rootTree.GetEntries()
print(f" ==> Start flatteing {nEntries} Events")

##################
### Flattening  ##
##################

#will hold the nanoAOD values -> can handle jagged arrays!!
nf = NanoFrame(fileName, ) 

#dictionary which holds everything
nfDir = {}
#dictionary which holds selected candidates
cands = {}
#dictionary which holds final candidate
finalCand = {}

#branches and names
branches = rootTree.GetListOfBranches()
names    = [branch.GetName() for branch in branches]
# remove the last few branches (f.e. fixedGridRhoFastjetAll) 
# Do we need them ? otherwise I have to adapt the shape as below
# Where do these Muon branches come from?? I have to find and remove
names = [name for name in names if not (name.startswith("fixed") or name == "nMuon" or name.startswith("Muon")) ] 

#Bs mass
bsMass_ = 5.36688

for name in names:

  #Define all your branches
  nfDir[name] = nf[name]
  cands[name] = [] #lets take lists because its less painful 
  finalCand[name] = []


#adapt these shapes to the nr of candidates in every event
nfDir["run"]             = np.array([np.array([run]  *len(nfDir["mu_pt"][ev])) for ev,run  in enumerate(nfDir["run"]) ],             dtype=object)
nfDir["luminosityBlock"] = np.array([np.array([lumi] *len(nfDir["mu_pt"][ev])) for ev,lumi in enumerate(nfDir["luminosityBlock"]) ], dtype=object)
nfDir["event"]           = np.array([np.array([evt]  *len(nfDir["mu_pt"][ev])) for ev,evt  in enumerate(nfDir["event"]) ],           dtype=object)
nfDir["n"]               = np.array([np.array([n]    *len(nfDir["mu_pt"][ev])) for ev,n    in enumerate(nfDir["n"]) ],               dtype=object)

#now we loop over the events and select the good candidates
for ev in range(nEntries):

  # for every event build a flag based on the selections
  # the flag will tell which candidates to drop *per event*, flag = np.array([True,False,True,...])

  # To discuss: what else?
  flag = ((nfDir["kk_charge"][ev]         < 0) * 
          (nfDir["gen_match_success"][ev] > 0) *
          (nfDir["dsMu_mass"][ev]         < bsMass_ ))
          # (nfDir["pi_mu_charge"] > 0 )

  for name in names:
    cands[name].append(nfDir[name][ev][flag]) 
    
  
  # decide to pick the candidate with highest dsMu pt -> discuss

  #check if not empty (if one var is empty, all of them are)
  if (cands["dsMu_pt"][ev].size > 0 ):
    # get max pt
    maxPtFlag = np.argmax(cands["dsMu_pt"][ev])
    # select the final candidate
    for name in names:
      finalCand[name].append(cands[name][ev][maxPtFlag])

"""
# Create an RDataFrame with an index column
rdf = ROOT.RDataFrame(nEntries)
# Define branches for each float in the container
for i, value in enumerate(floats):
    branch_name = f"float_{i}"
    rdf = rdf.Define(branch_name, f"double({value})")

# Write the RDataFrame to a ROOT tree in the output file
rdf.Snapshot("tree", "output_file.root")
"""
