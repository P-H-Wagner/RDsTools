import ROOT
import argparse
import os

from helper import colors, labels, bsMass
from histModels import models

# parsing
#parser = argparse.ArgumentParser()
#parser.add_argument('filename')
#args = parser.parse_args()

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

#########################################
## INPUT                               ##
#########################################

sig = "22_03_2024_17_34_45" #old signals 
hb  = "22_03_2024_17_34_53" #old inclusive

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
## CREATE RDF FROM TREE                ##
#########################################

def getRdf(dateTime):

  #access the flat ntuples
  files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test

  #chain them
  chain = ROOT.TChain("tree")
  chain.Add(files)

  #create rdf from tree
  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)

# Create rdf from tree
print(f" ===> Start creating RDataFrames")
chainSig, rdfSig = getRdf(sig)
chainHb,  rdfHb  = getRdf(hb)

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

  if "mutau" in key:    
    labels.append(r"#mu - signals");   labels.append(r"#tau - signals")
    colors.append(ROOT.kAzure);         colors.append(ROOT.kPink)
  
  elif "dsdsstar" in key: 
    labels.append(r"D_{s} - signals"); labels.append(r"D*_{s} - signals")
    colors.append(ROOT.kMagenta - 4); colors.append(ROOT.kGreen + 3)

  elif "allsignals" in key:      
    labels.append(r"D_{s} #mu"); labels.append(r"D_{s} #tau"); labels.append(r"D*_{s} #mu");   labels.append(r"D*_{s} #tau")
    colors.append(ROOT.kAzure);  colors.append(ROOT.kPink);    colors.append(ROOT.kAzure + 8); colors.append(ROOT.kPink + 5)

  if "hb" in key:
    labels.append("Hb")
    colors.append(ROOT.kGray + 2)

  return (colors,labels)

#########################################
## CREATE LEGEND                       ##
#########################################

def getLegend(x0 = 0.5, y0 = 0.8, x1 = 0.9, y1 = 0.9, nColumns = 1):

  "Default setting is x0,y0,x1,y1 as given with 1 column"

  legend = ROOT.TLegend(x0,y0,x1,y1)
  legend.SetTextSize(0.04)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetNColumns(nColumns) 
 
  return legend

#########################################
## GET YMAX FOR PLOTTING               ##
#########################################

def getYMax(histos, norm = True):

  maxis = []
  for hist in histos:
    if norm: hist.Scale(1/hist.Integral())
    maxis.append(hist.GetMaximum())
 
  return max(maxis)

#########################################
## TOOLS FOR WEIGHTED RECO METHOD      ##
#########################################


def getReco1Scale(selection, name,chain = chainSig):

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
## GINI COEFFICIENT                    ## 
#########################################

def getGini(histos,key):

  if key == "mutau": label = r"Gini #mu vs. #tau"
  else:              label = "D_{s} vs. D*_{s}" 

  areaTot, area1, area2 = 0,0,0

  for i in range(histos[0].GetNbinsX()):
   
    entry1 = histos[0].GetBinContent(i+1)
    entry2 = histos[1].GetBinContent(i+1)
    areaTot += entry1 + entry2
 
    if (entry1 >= entry2):
      area1 += entry1 - entry2 
    else:
      area2 += entry2 - entry1 

  giniCoeff = round((area1 + area2) / areaTot,2)

  giniText = ROOT.TPaveText(0.2, 0.8, 0.5, 0.9, 'nbNDC')
  giniText.AddText(f'{label} = {giniCoeff}')
  giniText.SetTextFont(42)
  giniText.SetTextSize(0.04)
  giniText.SetFillColor(0)
  giniText.SetFillStyle(0)
  giniText.SetBorderSize(0)

  return giniText

#########################################
## CREATE DEFAULT HISTOS               ##
#########################################

def createHistos(selection,rdf):

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

    try:    legend.AddEntry(hist,            labels[i],"l")
    except: legend.AddEntry(hist.GetValue(), labels[i],"l")
    legend.Draw("SAME")
 
    #get gini coefficietna as separation power
    giniText = getGini(histos[0:2], key)
    giniText.Draw("EP SAME")
  
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

def callSignalDistinction(histoType1, histoType2, selection):
  
  for var in models.keys():
  
    #mu tau only
    signalDistinction([histoType1[var],  histoType2[var]],                        var, "mutau")
    #mu tau and hb
    signalDistinction([histoType1[var],  histoType2[var], histosBasicHb[var]],    var, "mutauhb")
    #weighted reco method
  
    if "reco_1" in var:
  
      name  = var[:-7]
      reco1 = var
      reco2 = name + "_reco_2" 
  
      weightedMu  = getWeightedReco({reco1: histoType1[reco1],      reco2: histoType1[reco2]},      name, selection)
      weightedTau = getWeightedReco({reco1: histoType2[reco1],      reco2: histoType2[reco2]},      name, selection)
      weightedHb  = getWeightedReco({reco1: histosBasicHb[reco1],   reco2: histosBasicHb[reco2]},   name, selection)
  
      signalDistinction([weightedMu, weightedTau],             "reco_weighted", "mutau")
      signalDistinction([weightedMu, weightedTau, weightedHb], "reco_weighted", "mutauhb")
  


#########################################
## MAIN                                ##
#########################################

# Create histos for different selections
print(f" ===> Start creating histograms")

#Signal
histosBasic     = createHistos(selBasic,  rdfSig)
histosMu        = createHistos(selMu,     rdfSig)
histosTau       = createHistos(selTau,    rdfSig)
histosDs        = createHistos(selDs,     rdfSig)
histosDsStar    = createHistos(selDsStar, rdfSig)

#Hb
histosBasicHb     = createHistos(selBasic,  rdfHb)


print(f" ===> Basic, Mu, Tau, Ds, DsStar DONE")

#histosDsMu      = createHistos(selDsMu)
#histosDsTau     = createHistos(selDsTau)
#histosDsStarMu  = createHistos(selDsStarMu)
#histosDsStarTau = createHistos(selDsStarTau)
print(f" ===> Single Signals DONE")


# Plot Mu vs Tau
callSignalDistinction(histosMu, histosTau,    selBasic)
# Plot Ds vs Ds star
callSignalDistinction(histosDs, histosDsStar, selBasic)

#Method comparison

#prefix = ["q2","m2_miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt"]
prefix = ["q2","m2_miss"]
"""
for name in prefix:
  print(f" ===> Compare Bs reco methods for {name}")
  # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...
  histos = {var: histosDsMu[var] for var in models.keys() if name in var}
  methodDistinction(histos,name,selDsMu)

"""

