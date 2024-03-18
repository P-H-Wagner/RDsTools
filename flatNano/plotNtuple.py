import ROOT
import argparse
from helper import colors, labels, bsMass
from histModels import models

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# disable title and stats
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

#input files
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/*"
#files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/all_signals_flatChunk_0.root" #test

#chain them
chain = ROOT.TChain("tree")
chain.Add(files)

#create rdf from tree
rdf = ROOT.RDataFrame(chain)

#get branche names
names = [branch.GetName() for branch in chain.GetListOfBranches()]

#selections
selBasic        = "(bs_pt_reco_1 == bs_pt_reco_1) && (dsMu_mass < bsMass)"  # bs cut and reco are not nan

selDsMu         = "(sig == 0)"                                        # select Ds  Mu 
selDsTau        = "(sig == 1)"                                        # select Ds  Tau 
selDsStarMu     = "(sig == 5)"                                        # select Ds* Mu
selDsStarTau    = "(sig == 6)"                                        # select Ds* Tau

selMu           = "(sig == 0 || sig == 5)"                                  # select mu  signals
selTau          = "(sig == 1 || sig == 6)"                                  # select tau signals
selDs           = "(sig == 0 || sig == 1)"                                  # select Ds signals
selDsStar       = "(sig == 5 || sig == 6)"                                  # select Ds star signals

#Entries
nTot      = chain.GetEntries()
nBasic    = chain.GetEntries(selBasic)
nMu       = chain.GetEntries(selMu)
nTau      = chain.GetEntries(selTau)
nDs       = chain.GetEntries(selDs)
nDsStar   = chain.GetEntries(selDsStar)

def createHistos(selection):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"

  histos = {}
   
    for var,model in models.items(): 
      histos[var] = rdf.Filter(selection).Histo1D(model)

  return histos

createHistos(selBasic)
