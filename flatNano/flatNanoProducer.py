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

import argparse
parser = argparse.ArgumentParser(
                     prog='ProgramName',
                     description='What the program does',
                     epilog='Text at the bottom of help')
 
parser.add_argument('filename')
parser.add_argument('channel')

args = parser.parse_args()

#Single File
#fileName = "/work/pahwagne/test/multi_instances.root" 
#load it into root file to get all branch names
#rootFile = ROOT.TFile.Open(fileName)
#rootTree = rootFile.Get("Events")

#before we put the bs reconstruction methods into functions
#inp = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/11_02_2024_21_56_13/"
#and after


inp = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/{args.filename}/"

out =  f"/scratch/pahwagne/flatNano2"
#out = f"/work/pahwagne/test/flatNano/{args.channel}/" #used for local testy only

files = os.listdir(inp)
files = files[0:1]

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
  toSave = {}
  
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
    toSave[name] = []
  
  
  #adapt these shapes to the nr of candidates in every event
  nfDir["run"]             = np.array([np.array([run]  *len(nfDir[names[-1]][ev])) for ev,run  in enumerate(nfDir["run"]) ],             dtype=object)
  nfDir["luminosityBlock"] = np.array([np.array([lumi] *len(nfDir[names[-1]][ev])) for ev,lumi in enumerate(nfDir["luminosityBlock"]) ], dtype=object)
  nfDir["event"]           = np.array([np.array([evt]  *len(nfDir[names[-1]][ev])) for ev,evt  in enumerate(nfDir["event"]) ],           dtype=object)
  nfDir["n"]               = np.array([np.array([n]    *len(nfDir[names[-1]][ev])) for ev,n    in enumerate(nfDir["n"]) ],               dtype=object)

  
  #now we loop over the events and select the good candidates
  for ev in range(1): #nEntries):
  
    #get the nr of candidates in the event
    nCand = len(nfDir[names[-1]][ev])

    #collect all candidates per event in a list
    allCandidates = []

    #test = [] # for debugging

    for iCand in range(nCand):
      # define empty list for every candidate
      dummy = []
      #dummy2 = [] # for debugging
      for i,name in enumerate(names):
        dummy.append(nfDir[name][ev][iCand]) 
        # keep important indices for sorting
        if (name == "mu_charge"): iMuCharge = i    
        if (name == "pi_charge"): iPiCharge = i    
        if (name == "k1_charge"): iK1Charge = i    
        if (name == "k2_charge"): iK2Charge = i    
        if (name == "dsMu_mass"): iDsMuMass = i    
        if (name == "dsMu_pt")  : iDsMuPt   = i    
     
      allCandidates.append(dummy)

      #dummy2.append(dummy[iMuCharge]) #for debugging
      #dummy2.append(dummy[iPiCharge])
      #dummy2.append(dummy[iK1Charge])
      #dummy2.append(dummy[iK2Charge])
      #dummy2.append(dummy[iDsMuMass])
      #dummy2.append(dummy[iDsMuPt])
      #test.append(dummy2) 
     
    # sort the candidates according to:
    # 1. we want opposite pi mu charge (priority)
    # 2. we want opposite kk charge (2nd priority)
    # 3. we want the dsMu mass to be physical (lower than Bs mass)
    # 4. we sort after dsMu pt

    # The logic goes as follows: if we have two candidates, one has the correct dsMu charge
    # assignment but low dsMu pt, we still keep this event over wrong charge candidates with
    # high pt, according to the priority list above. We do not throw away events! Worst case
    # scenario is that we have only wrong charge candidates and take the one with highest pt
   
 
    allCandidates.sort(key = lambda x : ( (x[iMuCharge]*x[iPiCharge] < 0.,x[iK1Charge]*x[iK2Charge] < 0. , x[iDsMuMass] < bsMass_, x[iDsMuPt]) ), reverse = True) 

    # the best candidate is the first one
    bestCandidate = allCandidates[0]

    # save the best candidate
    for i,name in enumerate(names):
        toSave[name].append(bestCandidate[i])
  
  # events which survived (take one column) should be the same -> checked
  nFinal = len(toSave[names[-1]])
  
  # output file
  output_file = ROOT.TFile( out + f"/{args.channel}_flatNano_chunk_{chunkNr}.root", "RECREATE")
  
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
      containers[name][0] = toSave[name][i]
  
    # very important to call Fill() outside the name loop!
    tree.Fill()
  
  # Write the TTree to the ROOT file
  output_file.Write()
  # Close the ROOT file
  output_file.Close()


