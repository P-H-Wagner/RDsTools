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
import timeit


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

#out =  f"/scratch/pahwagne/flatNano2"
out = f"/work/pahwagne/test/flatNano/{args.channel}/" #used for local testy only

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
  nEntries = 100
  print(f" ==> Start flatteing {nEntries} Events of chunk {chunkNr}")
  
  ##################
  ### Flattening  ##
  ##################
  
  #will hold the nanoAOD values -> can handle jagged arrays!!
  nf = NanoFrame(fileName, ) 
  
  #dictionary which holds everything
  nfDir = {}
  #list which holds number of candidates in every event
  nCands = []
  #dictionary which holds selected candidates
  cands = {}
  #dictionary which holds final candidate
  toSave = {}
  
  #branches and names
  branches = rootTree.GetListOfBranches()
  names    = [branch.GetName() for branch in branches]
  # remove the last few branches (f.e. fixedGridRhoFastjetAll) 
  # Do we need them ? otherwise I have to adapt the shape as below
  names = [name for name in names if not (name.startswith("fixed") or name == "nmuon" or name.startswith("muon")) ] 
  
  #Bs mass
  bsMass_ = 5.36688

  print(f"We take over the shape from: {names[-1]}")

  for i,name in enumerate(names):

    #Define all your branches
    nfDir[name] = nf[name]
    cands[name] = [] #lets take lists because its less painful 
    toSave[name] = []

    # keep important indices for sorting
    if (name == "mu_charge"): iMuCharge = i    
    if (name == "pi_charge"): iPiCharge = i    
    if (name == "k1_charge"): iK1Charge = i    
    if (name == "k2_charge"): iK2Charge = i    
    if (name == "dsMu_mass"): iDsMuMass = i    
    if (name == "dsMu_pt")  : iDsMuPt   = i    

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
  nfDir["run"]             = np.array([np.array([run]  *len(nfDir[names[-1]][ev])) for ev,run  in enumerate(nfDir["run"]) ],             dtype=object)
  nfDir["luminosityBlock"] = np.array([np.array([lumi] *len(nfDir[names[-1]][ev])) for ev,lumi in enumerate(nfDir["luminosityBlock"]) ], dtype=object)
  nfDir["event"]           = np.array([np.array([evt]  *len(nfDir[names[-1]][ev])) for ev,evt  in enumerate(nfDir["event"]) ],           dtype=object)
  nfDir["n"]               = np.array([np.array([n]    *len(nfDir[names[-1]][ev])) for ev,n    in enumerate(nfDir["n"]) ],               dtype=object)
  try: nfDir["ngen"]       = np.array([np.array([n]    *len(nfDir[names[-1]][ev])) for ev,n    in enumerate(nfDir["ngen"]) ],            dtype=object)
  except: pass

  #start timing
  start = timeit.timeit()
   

  # Lets rewrite the nfDir in the following shape and sort at the same time to avoid an extra loop:
  #          [        [ [this is one cand with all its variables pt, mass, m2miss, .... ]        #and we dump all candidates in an array]  sort this!                                # for each event]          
  allCands = [ sorted([ [nfDir[name][ev][nCand] for name in names                       ] for nCand in range(len(nfDir[names[-1]][ev])) ], key = mySorting, reverse = True) for ev in range(nEntries)]

  # Lets sort them for each event (each element in allCands, as this is the outermost index (see above) and pick the best (index 0)
  bestCands = np.array([ listOfCands[0] for listOfCands in allCands ]) 
 
  ## For saving, we have to slice after the name rather than event
  for i,name in enumerate(names):
    toSave[name] = bestCands[:,i]

  outfile = uproot.recreate("example.root")
  outfile["tree"] = pd.DataFrame(toSave)

  """ 
  #loop over the variable names
  for name in names:
    #collect all candidates per event in a list (a list per event!)
    allCands[name]       = [ list(nf[name][ev])  for ev in range(nEntries)] 
    #sort the candidates for each event (sorted list per event!)
  

  allCandsSorted[name] = [ candList.sort(key = lambda x : ( (x[iMuCharge]*x[iPiCharge] < 0.,x[iK1Charge]*x[iK2Charge] < 0. , x[iDsMuMass] < bsMass_, x[iDsMuPt]) ), reverse = True) for candList in allCandidates ]
    #pick the best candidate per event (one nr per event!)
    bestCandidates      = [ candList[0] for candList in allCandidatesSorted]
    #save the best candidates
    toSave[name]        = bestCandidates 
     
    # sort the candidates according to:
    # 1. we want opposite pi mu charge (priority)
    # 2. we want opposite kk charge (2nd priority)
    # 3. we want the dsMu mass to be physical (lower than Bs mass)
    # 4. we sort after dsMu pt

    # The logic goes as follows: if we have two candidates, one has the correct dsMu charge
    # assignment but low dsMu pt, we still keep this event over wrong charge candidates with
    # high pt, according to the priority list above. We do not throw away events! Worst case
    # scenario is that we have only wrong charge candidates and take the one with highest pt
  stop = timeit.timeit() 

  print(f"used time: {stop - start}")
  # events which survived (take one column) should be the same -> checked
  nFinal = len(toSave[names[-1]])
  """
  """ 
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
  """

