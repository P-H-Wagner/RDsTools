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

sig  = "09_04_2024_18_02_43"  #"22_03_2024_17_34_45" #old signals 
hb   = "09_04_2024_18_44_36"  #"22_03_2024_17_34_53" #old inclusive
data = "" #data -> apply SB method
#########################################
## SELECTIONS                          ##
#########################################

#selections
selBasic        = f"(disc_is_negative == 0) && (dsMu_mass < {bsMass})"                    # bs cut and disc not neg
selBasic        = f"dsMu_mass < {bsMass}"
selBasicHb      = selBasic + "&& (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 5) && (gen_sig != 6)"    # exlude signals from hb

selDsMu         = selBasic + "&& (gen_sig == 0)"                                              # select Ds  Mu 
selDsTau        = selBasic + "&& (gen_sig == 1)"                                              # select Ds  Tau 
selDsStarMu     = selBasic + "&& (gen_sig == 5)"                                              # select Ds* Mu
selDsStarTau    = selBasic + "&& (gen_sig == 6)"                                              # select Ds* Tau

selMu           = selBasic + "&& (gen_sig == 0 || gen_sig == 5)"                                  # select mu  signals
selTau          = selBasic + "&& (gen_sig == 1 || gen_sig == 6)"                                  # select tau signals
selDs           = selBasic + "&& (gen_sig == 0 || gen_sig == 1)"                                  # select Ds signals
selDsStar       = selBasic + "&& (gen_sig == 5 || gen_sig == 6)"                                  # select Ds star signals

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

  if   "gen"      in var: color,label = ROOT.kMagenta,    "Gen"
  elif "coll"     in var: color,label = ROOT.kCyan,       "Coll."
  elif "lhcb_alt" in var: color,label = ROOT.kViolet,     "LHCb xyz"
  elif "lhcb"     in var: color,label = ROOT.kBlue,       "LHCb z"
  elif "reco_1"   in var: color,label = ROOT.kGreen,      "Reco 1"
  elif "reco_2"   in var: color,label = ROOT.kOrange,     "Reco 2"
  elif "weighted" in var: color,label = ROOT.kYellow + 2, "Reco Weighted"

  else:                   color,label = ROOT.kBlack,   "Reco"

  return (color,label)

def getHbColorAndLabel():
  return (ROOT.kGray + 2, "Hb")


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
    colors.append(ROOT.kAzure);  colors.append(ROOT.kPink);    colors.append(ROOT.kMagenta - 4); colors.append(ROOT.kGreen + 3)

  if "hb" in key:
    hbColor,hbLabel = getHbColorAndLabel()
    labels.append(hbLabel)
    colors.append(hbColor)

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
 
  if not isinstance(histos, list): histos = [histos]

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
  nReco1 = chain.GetEntries(selection + f"&& (abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) ")
    
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

  if "mutau" in key: label    = r"Gini #mu vs. #tau"
  if "dsdsstar" in key: label = r"Gini D_{s} vs. D*_{s}" 
  if "sighb" in key:  label   = r"Gini Sig. vs. Hb"
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

def simpleDraw(histos, var, norm = True, key="simple"):

  #only makes directory if not already existing
  os.system(f"mkdir -p ./{key}/")

  yMax = getYMax(histos)

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(nColumns = 2)
  canvas.cd()

  # get default labels
  _,label = getColorAndLabel(var)

  if norm: histos[0].Scale(1/histos[0].Integral())

  histos[0].SetMaximum(1.3*yMax)
  histos[0].Draw("HIST")
  try: legend.AddEntry(histos[0].GetValue(),label, "l")
  except: legend.AddEntry(histos[0],label, "l")

  if "hb" in key: 
      
     giniText = getGini(histos, key = "sighb")
     giniText.Draw("EP SAME")

     if norm: histos[1].Scale(1/histos[1].Integral())
     hbColor,hbLabel = getHbColorAndLabel()
     histos[1].SetLineColor(hbColor)
     histos[1].Draw("SAME HIST")
     try: legend.AddEntry(histos[1].GetValue(),hbLabel, "l")
     except: legend.AddEntry(histos[1],hbLabel, "l")
  legend.Draw("SAME")

  canvas.SaveAs(f"./{key}/{key}_{var}.png")  

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
    if len(histos) < 4:
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

def callSignalDistinction(histos, selections, key):
  
  for var in models.keys():
 
    if len(histos) == 3:
      # we compare 2 types only

      #type1 type2 only
      signalDistinction([histos[0][var],  histos[1][var]],                        var, key)
      #inclduing hb
      signalDistinction([histos[0][var],  histos[1][var], histos[2][var]],    var, key+"hb")
      #weighted reco method
    
      if "reco_1" in var:
    
        name  = var[:-7]
        reco1 = var
        reco2 = name + "_reco_2" 
    
        weightedMu  = getWeightedReco({reco1: histos[0][reco1],      reco2: histos[0][reco2]},      name, selections[0])
        weightedTau = getWeightedReco({reco1: histos[1][reco1],      reco2: histos[1][reco2]},      name, selections[0])
        weightedHb  = getWeightedReco({reco1: histos[2][reco1],      reco2: histos[2][reco2]},      name, selections[1])
    
        signalDistinction([weightedMu, weightedTau],             name + "_reco_weighted", key)
        signalDistinction([weightedMu, weightedTau, weightedHb], name + "_reco_weighted", key+"hb")

    else:
      #all signals
      signalDistinction([histos[0][var],  histos[1][var], histos[2][var], histos[3][var], histos[4][var]],    var, key+"hb")

      if "reco_1" in var:
    
        name  = var[:-7]
        reco1 = var
        reco2 = name + "_reco_2" 
    
        weightedDsMu       = getWeightedReco({reco1: histos[0][reco1],      reco2: histos[0][reco2]},      name, selections[0])
        weightedDsTau      = getWeightedReco({reco1: histos[1][reco1],      reco2: histos[1][reco2]},      name, selections[0])
        weightedDsStarMu   = getWeightedReco({reco1: histos[2][reco1],      reco2: histos[2][reco2]},      name, selections[0])
        weightedDsStarTau  = getWeightedReco({reco1: histos[3][reco1],      reco2: histos[3][reco2]},      name, selections[0])
        weightedHb         = getWeightedReco({reco1: histos[4][reco1],      reco2: histos[4][reco2]},      name, selections[1])


        signalDistinction([weightedDsMu, weightedDsTau, weightedDsStarMu, weightedDsStarTau, weightedHb],             name + "_reco_weighted", key+"hb")

  
def callSimpleDraw(histosSig,histosHb, selection,selectionHb):

  for var in models.keys():
 
    simpleDraw([histosSig[var], histosHb[var]], var, key = "sighb")

    if "reco_1" in var:
  
      name  = var[:-7]
      reco1 = var
      reco2 = name + "_reco_2" 
  
      weightedMu  = getWeightedReco({reco1: histosSig[reco1],       reco2: histosSig[reco2]},       name, selection)
      weightedHb  = getWeightedReco({reco1: histosHb[reco1],        reco2: histosHb[reco2]},        name, selectionHb)
  
      simpleDraw([weightedMu, weightedHb], name + "_reco_weighted", key = "sighb")



#########################################
## MAIN                                ##
#########################################

# Create histos for different selections
print(f" ===> Start creating histograms")

#Signal
#histosBasic     = createHistos(selBasic,  rdfSig)
#histosMu        = createHistos(selMu,     rdfSig)
#histosTau       = createHistos(selTau,    rdfSig)
#histosDs        = createHistos(selDs,     rdfSig)
#histosDsStar    = createHistos(selDsStar, rdfSig)

#Hb
histosBasicHb     = createHistos(selBasicHb,  rdfHb)


print(f" ===> Basic, Mu, Tau, Ds, DsStar DONE")

histosDsMu      = createHistos(selDsMu,      rdfSig)
histosDsTau     = createHistos(selDsTau,     rdfSig)
histosDsStarMu  = createHistos(selDsStarMu,  rdfSig)
histosDsStarTau = createHistos(selDsStarTau, rdfSig)
print(f" ===> Single Signals DONE")

# Simple Plots 
#callSimpleDraw(histosBasic,histosBasicHb, selBasic,selBasicHb)

# Plot Mu vs Tau
#callSignalDistinction( [histosMu, histosTau, histosBasicHb],    [selBasic,selBasicHb], "mutau")
# Plot Ds vs Ds star
#callSignalDistinction( [histosDs, histosDsStar, histosBasicHb], [selBasic,selBasicHb], "dsdsstar")
# Plot all signals
callSignalDistinction( [histosDsMu, histosDsTau, histosDsStarMu, histosDsStarTau, histosBasicHb], [selBasic,selBasicHb], "allsignals")
#Method comparison

#prefix = ["q2","m2_miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt"]
"""
for name in prefix:
  print(f" ===> Compare Bs reco methods for {name}")
  # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...
  histos = {var: histosDsMu[var] for var in models.keys() if name in var}
  methodDistinction(histos,name,selDsMu)

"""

