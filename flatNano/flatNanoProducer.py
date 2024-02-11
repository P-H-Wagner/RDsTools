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
from array import array
import os
import sys

#Single File
#fileName = "/work/pahwagne/test/multi_instances.root" 
#load it into root file to get all branch names
#rootFile = ROOT.TFile.Open(fileName)
#rootTree = rootFile.Get("Events")


inp = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/11_02_2024_21_56_13/"
out =  "/scratch/pahwagne/flatNano"  
#out = "/work/pahwagne/test/flatNano" #used for local testy only

files = os.listdir(inp)
#files = files[0:3]

for f in files:
  fileName = inp + f
  chunkNr = fileName[-7:-5]

  #load it into root file to get all branch names
  rootFile = ROOT.TFile.Open(fileName)
  rootTree = rootFile.Get("Events")

  # nr of events
  nEntries = rootTree.GetEntries()
  print(f" ==> Start flatteing {nEntries} Events of chunk {chunkNr}")
  
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
  nfDir["run"]             = np.array([np.array([run]  *len(nfDir[names[-1]][ev])) for ev,run  in enumerate(nfDir["run"]) ],             dtype=object)
  nfDir["luminosityBlock"] = np.array([np.array([lumi] *len(nfDir[names[-1]][ev])) for ev,lumi in enumerate(nfDir["luminosityBlock"]) ], dtype=object)
  nfDir["event"]           = np.array([np.array([evt]  *len(nfDir[names[-1]][ev])) for ev,evt  in enumerate(nfDir["event"]) ],           dtype=object)
  nfDir["n"]               = np.array([np.array([n]    *len(nfDir[names[-1]][ev])) for ev,n    in enumerate(nfDir["n"]) ],               dtype=object)

  
  #now we loop over the events and select the good candidates
  for ev in range(nEntries):
  
    # for every event build a flag based on the selections
    # the flag will tell which candidates to drop *per event*, flag = np.array([True,False,True,...])
 
    # To discuss: what else?
    flag = ((nfDir["kk_charge"][ev]        == 0) * #future samples: k_k_charge < 0 
            (nfDir["gen_match_success"][ev] > 0) *
            (nfDir["dsMu_mass"][ev]         < bsMass_ ))
            # (nfDir["pi_mu_charge"] > 0 )
 
    print(flag) 
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
  
  # events which survived (take one column)
  nFinal = len(finalCand[names[-1]])
  
  # output file
  output_file = ROOT.TFile( out + f"/flatNano_chunk_{chunkNr}.root", "RECREATE")
  
  # crate tree
  tree = ROOT.TTree("tree", "tree")
  
  # container
  containers = {}
  
  for name in names:
    # container that will hold the values to fill in the tree 
    containers[name] = array('f', [0])
  
    # create branch
    tree.Branch(name, containers[name], f"{name}/F")  # "F" specifies that it's a float
  
  # fill the branch 
  for i in range(nFinal):
    for name in names:
      containers[name][0] = finalCand[name][i]
  
    # very important to call Fill() outside the name loop!
    tree.Fill()
  
  # Write the TTree to the ROOT file
  output_file.Write()
  # Close the ROOT file
  output_file.Close()


