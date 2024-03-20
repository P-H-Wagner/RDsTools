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


########################################
# TODO: include bkg !


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

selDsMu         = selBasic + "&& (sig == 0)"                                              # select Ds  Mu 
selDsTau        = selBasic + "&& (sig == 1)"                                              # select Ds  Tau 
selDsStarMu     = selBasic + "&& (sig == 5)"                                              # select Ds* Mu
selDsStarTau    = selBasic + "&& (sig == 6)"                                              # select Ds* Tau

selMu           = selBasic + "&& (sig == 0 || sig == 5)"                                  # select mu  signals
selTau          = selBasic + "&& (sig == 1 || sig == 6)"                                  # select tau signals
selDs           = selBasic + "&& (sig == 0 || sig == 1)"                                  # select Ds signals
selDsStar       = selBasic + "&& (sig == 5 || sig == 6)"                                  # select Ds star signals

#########################################
## COLOR SCHEME                        ##
#########################################

def getColorAndLabel(var):

  "Colors for simple plot of signgle histograms, default" 

  if   "gen"      in var: color,label = ROOT.kMagenta, "Gen"
  elif "coll"     in var: color,label = ROOT.kBlack,   "Coll."
  elif "lhcb_alt" in var: color,label = ROOT.kBlue,    "LHCb xyz"
  elif "lhcb"     in var: color,label = ROOT.kRed,     "LHCb z"
  elif "reco_1"   in var: color,label = ROOT.kGreen,   "Reco 1"
  elif "reco_2"   in var: color,label = ROOT.kOrange,  "Reco 2"

  else:                   color,label = ROOT.kGray+2,  "Reco"

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
  legend.SetNColumns(nColumns) 
 
  return legend

def getYMax(histos, norm = True):

  maxis = []
  for hist in histos:
    if norm: hist.Scale(1/hist.Integral())
    maxis.append(hist.GetMaximum())
 
  return max(maxis)

def getReco1Scale(selection, name):

  ntot   = chain.GetEntries(selection)
  nReco1 = chain.GetEntries(selection + f"&& (abs({name}_reco_1 - {name}_gen) < abs({name}_reco_2 - {name}_gen)) ")
    
  return nReco1 / ntot 

def getWeightedReco(histos, name, selection, norm = True):

  #get reco scale
  reco1Scale = getReco1Scale(selection,name)

  #clone recos
  reco1Copy = histos[name + "_reco_1"].Clone()
  reco2Copy = histos[name + "_reco_2"].Clone()

  #normalize
  if norm: 
    reco1Copy.Scale(1/histos[name + "_reco_1"].Integral())
    reco2Copy.Scale(1/histos[name + "_reco_2"].Integral())

    reco1Copy.GetYaxis().SetTitle("a.u.")

  #scale and add them
  reco1Copy.Scale(reco1Scale) 
  reco2Copy.Scale(1 - reco1Scale) 
  reco1Copy.Add(reco2Copy)
    
  return reco1Copy 


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

   print(f"      Selection: {selection} DONE")

  return histos

#########################################
## DRAWING FUNCTIONS                   ##
#########################################

########################################
# Simple, plots all histo models

def simpleDraw(histos):

  #only makes directory if not already existing
  os.system(f"mkdir -p ./simpleDraw/")

  for var in  histos.keys():
   
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    legend = getLegend()
    canvas.cd()

    # get default labels
    _,label = getColorAndLabel(var)

    histos[var].Draw("HIST")

    legend.AddEntry(histos[var].GetValue(),label, "l")
    legend.Draw("SAME")

    canvas.SaveAs(f"./simpleDraw/{var}.png")  

########################################
# Plots signal distinction according to key 

def signalDistinction(histos,var,key,norm = True):

  "Histos is a LIST!"

  #only makes directory if not already existing
  os.system(f"mkdir -p ./{key}/")

  # get new colors and labels for distinction
  colors, labels = getColorAndLabelSignalDistinction(key) 

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(nColumns = 2)
  canvas.cd()

  # get maximum
  yMax = getYMax(histos,norm)

  for i, hist in enumerate(histos):
    
    hist.SetLineColor(colors[i]) #adapt color
    if norm: 
      hist.Scale(1/hist.Integral())
      hist.GetYaxis().SetTitle("a.u.")
    hist.SetMaximum(1.3*yMax)

    # draw 
    if i == 0: hist.Draw("HIST")
    else:      hist.Draw("HIST SAME")

    legend.AddEntry(hist.GetValue(),labels[i],"l")
    legend.Draw("SAME")
 
  
  canvas.SaveAs(f"./{key}/{key}_{var}.png")  

########################################
# Plots bs reconstruction method distinction 

def methodDistinction(histos,name, selection,norm = True):

  "Histos is a DICTIONARY!"

  #only makes directory if not already existing
  os.system(f"mkdir -p ./method_dist/")

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(nColumns = 3)
  canvas.cd()

  yMax = getYMax(histos.values(),norm)

  #recreate the weighted combination
  recoAvailable = name + "_reco_1" in models.keys() #check that recos are available 
    
  if recoAvailable: 
    weightedReco = getWeightedReco(histos,name,selection,norm)
    weightedReco.Draw("HIST") 
    weightedReco.SetMaximum(1.3*yMax)

  for i, hist in enumerate(histos.values()):

    var = list(histos.keys())[i]

    if norm: 
      hist.Scale(1/hist.Integral())
      hist.GetYaxis().SetTitle("a.u.")

    hist.SetMaximum(1.3*yMax)    
    _,label = getColorAndLabel(var)

    # draw 
    if ((i == 0) and (not recoAvailable)): hist.Draw("HIST")
    else:                                  hist.Draw("HIST SAME")

    legend.AddEntry(hist.GetValue(),label,"l")
    legend.Draw("SAME")

  canvas.SaveAs(f"./method_dist/method_dist_{name}.png")  
  print(f"DONE")


#########################################
## MAIN                                ##
#########################################


# Create histos for different selections
print(f" ===> Start creating histograms")

histosBasic     = createHistos(selBasic)
histosMu        = createHistos(selMu)
histosTau       = createHistos(selTau)
histosDs        = createHistos(selDs)
histosDsStar    = createHistos(selDsStar)
print(f" ===> Basic, Mu, Tau, Ds, DsStar DONE")

histosDsMu      = createHistos(selDsMu)
histosDsTau     = createHistos(selDsTau)
histosDsStarMu  = createHistos(selDsStarMu)
histosDsStarTau = createHistos(selDsStarTau)
print(f" ===> Single Signals DONE")

#simpleDraw(histosTau, "tau_only")

# Plot Mu vs Tau
"""
for var in models.keys():
  signalDistinction([histosMu[var],histosTau[var]], var, "mutau")
# Plot Ds vs DsStar
for var in models.keys():
  signalDistinction([histosDs[var],histosDsStar[var]], "dsdsstar")

# Plot all signals (Ds mu, Ds tau, Ds* mu Ds* tau)
for var in models.keys():
  signalDistinction([histosDsMu[var],histosDsTau[var],histosDsStarMu[var],histosDsStarTau[var]], "allsignals")
"""
#Method comparison

#prefix = ["q2","m2_miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt"]
prefix = ["q2","m2_miss"]

for name in prefix:
  print(f" ===> Compare Bs reco methods for {name}")
  # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...
  histos = {var: histosDsMu[var] for var in models.keys() if name in var}
  methodDistinction(histos,name,selDsMu)



