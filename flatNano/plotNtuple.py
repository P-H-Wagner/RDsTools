import ROOT
import argparse
import os

from helper import colors, labels, bsMass
from histModels import models

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)


#########################################
## INPUT                               ##
#########################################

#input files
#files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/*"
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/all_signals_flatChunk_0.root" #test

#chain them
chain = ROOT.TChain("tree")
chain.Add(files)

#create rdf from tree
rdf = ROOT.RDataFrame(chain)

#get branche names
names = [branch.GetName() for branch in chain.GetListOfBranches()]

#########################################
## SELECTIONS                          ##
#########################################

#selections
selBasic        = f"(bs_pt_reco_1 == bs_pt_reco_1) && (dsMu_mass < {bsMass})"  # bs cut and reco are not nan

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

#########################################
## COLOR SCHEME                        ##
#########################################

def getColorAndLabel(var):

  "Colors for simple plot of signgle histograms, default" 

  if   "gen"      in var: color = ROOT.kMagenta; label = "Gen"
  elif "coll"     in var: color = ROOT.kBlack;   label = "Coll."
  elif "lhcb_alt" in var: color = ROOT.kBlue;    label = "LHCb xyz"
  elif "lhcb"     in var: color = ROOT.kRed;     label = "LHCb z"
  elif "reco_1"   in var: color = ROOT.kGreen;   label = "Reco 1"
  elif "reco_2"   in var: color = ROOT.kOrange;  label = "Reco 2"

  else: color = ROOT.kGray + 2; label = "Reco"

  return (color,label)

def getColorAndLabelSignalDistinction(key):

  "Adapt colors for multi histogram plot, f.e. mu vs tau"

  labels = []
  colors = []

  if   key == "mutau":    
    labels.append(r"#mu - signals");   labels.append(r"#tau - signals")
    colors.append(ROOT.kAzure);         colors.append(ROOT.kPink)
  
  elif key == "dsdsstar": 
    labels.append(r"D_{s} - signals"); labels.append(r"D*_{s} - signals")
    colors.append(ROOT.kMagenta - 4); labels.append(ROOT.kGreen + 3)

  elif key == "allsignals":      
    labels.append(r"D_{s} #mu"); labels.append(r"D_{s} #tau"); labels.append(r"D*_{s} #mu");   labels.append(r"D*_{s} #tau")
    colors.append(ROOT.kAzure);  colors.append(ROOT.kPink);    colors.append(ROOT.kAzure + 8); colors.append(ROOT.kPink + 5)

  return (colors,labels)

def getLegend(x0 = 0.5, y0 = 0.75, x1 = 0.9, y1 = 0.9, nColumns = 1):

  "Default setting is x0,y0,x1,y1 as given with 1 column"

  legend = ROOT.TLegend(x0,y0,x1,y1)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  
  return legend

def getReco1Scale(selection, name):

  ntot   = chain.GetEntries(selBasic)
  nReco1 = chain.GetEntries(selBasic + f"&& (abs({name}_reco_1 - {name}_gen) < abs({name}_reco_2 - {name}_gen)) ")
    
  return nReco1 / ntot 



#########################################
## CREATE DEFAULT HISTOS               ##
#########################################

def createHistos(selection):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"

  histos = {}
   
  for var,model in models.items(): 

      histos[var] = rdf.Filter(selection).Histo1D(model[0], var)
      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetMaximum())

      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(2)

  return histos

#########################################
## DRAWING FUNCTIONS                   ##
#########################################

def simpleDraw(histos, subfolder):

  #only makes directory if not already existing
  #os.system(f"mkdir -p simpleDraw/{subfolder}")

  for key in  histos.keys():
   
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    legend = getLegend()
    canvas.cd()

    # get default labels
    _,label = getColorAndLabel(key)

    histos[key].Draw("HIST")

    legend.AddEntry(histos[key],label, "l")
    legend.Draw("SAME")

    canvas.SaveAs(f"./simpleDraw/{key}.png")  

def signalDistinction(histos,name,key):

  # get new colors and labels for distinction
  colors, labels = getColorAndLabelSignalDistinction(key) 

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(nColumns = len(histos)%2)
  canvas.cd()

  for i, hist in enumerate(histos):
    
    # adapt color
    hist.SetLineColor(colors[i])
    hist.Scale(1/hist.Integral())

    # draw 
    if i == 0: hist.Draw("HIST")
    else:      hist.Draw("HIST SAME")

    #legend.AddEntry(hist,labels[i],"l")
    #legend.Draw("SAME")

  canvas.SaveAs(f"./mu_vs_tau/{name}_{key}.png")  

def methodDistinction(histos,name):

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(nColumns = len(histos)%2)
  canvas.cd()

  #recreate the weighted combination
  recoAvailable = name + "_reco_1" in models.keys() #check that recos are available 
    
  if recoAvailable:

    reco1Scale = getReco1Scale(name)

    reco1Copy = histos[name + "_reco_1"].Clone()
    reco1Copy.Scale(reco1Scale) 
    reco2Copy = histos[name + "_reco_2"].Clone()
    reco2Copy.Scale(1 - reco1Scale) 
    reco1Copy.Add(reco2Copy)

    reco1Copy.Draw("HIST") 
     

  for i, hist in enumerate(histos):
    
    # adapt color
    hist.SetLineColor(colors[i])

    # draw 
    if ((i == 0) and (not recoAvailable)): hist.Draw("HIST")
    else:                         hist.Draw("HIST SAME")

    #legend.AddEntry(hist,labels[i],"l")
    #legend.Draw("SAME")

  canvas.SaveAs(f"./method_distinction/method_{name}.png")  


#########################################
## MAIN                                ##
#########################################


# Create histos for different selections
histosBasic     = createHistos(selBasic)
histosMu        = createHistos(selMu)
histosTau       = createHistos(selTau)
histosDs        = createHistos(selDs)
histosDsStar    = createHistos(selDsStar)
histosDsMu      = createHistos(selDsMu)
histosDsTau     = createHistos(selDsTau)
histosDsStarMu  = createHistos(selDsStarMu)
histosDsStarTau = createHistos(selDsStarTau)

#simpleDraw(histosTau, "tau_only")

# Plot Mu vs Tau
for var in models.keys():
  signalDistinction([histosMu[var],histosTau[var]], var, "mutau")
"""
# Plot Ds vs DsStar
for var in models.keys():
  signalDistinction([histosDs[var],histosDsStar[var]], "dsdsstar")

# Plot all signals
for var in models.keys():
  signalDistinction([histosDsMu[var],histosDsTau[var],histosDsStarMu[var],histosDsStarTau[var]], "allsignals")

#Method comparison

prefix = ["q2","m2miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt"]

for name in prefix:
  # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...

  histos = [histosDsMu[var] if name in var for var in models.keys()]
  
  methodDistinction(histos,name)
"""



