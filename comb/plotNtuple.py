import ROOT
import argparse
import os
from sidebands import getSigma, getABCS
from helper import colors, labels, bsMass
from histModels import models
import numpy as np
from cms_style import CMS_lumi

# parsing
#parser = argparse.ArgumentParser()
#parser.add_argument('filename')
#args = parser.parse_args()

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()

#########################################
## INPUT                               ##
#########################################

sig  = "26_04_2024_16_28_22" #"10_04_2024_00_22_57"#"09_04_2024_18_02_43"  #"22_03_2024_17_34_45" #old signals 
hb   = "26_04_2024_16_28_58" #"10_04_2024_00_32_53"#"09_04_2024_18_44_36"  #"22_03_2024_17_34_53" #old inclusive
data = "26_04_2024_18_08_15" #data -> apply SB method

n    = 5 #how many sigma until the sb
dsMass_ = 1.96834
bsMass_ = 5.36688
#########################################
## SELECTIONS                          ##
#########################################


## BASIC

selBasic        = f" (dsMu_m < {bsMass}) & (k1_charge*k2_charge < 0) & (mu_charge*pi_charge <0) & (gen_match_success == 1))"
selBasicHb      = selBasic + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 5) && (gen_sig != 6) "    # exlude signals from hb

selDsMu         = selBasic + " && (gen_sig == 0)"                                              # select Ds  Mu 
selDsTau        = selBasic + " && (gen_sig == 1)"                                              # select Ds  Tau 
selDsStarMu     = selBasic + " && (gen_sig == 5)"                                              # select Ds* Mu
selDsStarTau    = selBasic + " && (gen_sig == 6)"                                              # select Ds* Tau

selMu           = selBasic + " && (gen_sig == 0 || gen_sig == 5)"                                  # select mu  signals
selTau          = selBasic + " && (gen_sig == 1 || gen_sig == 6)"                                  # select tau signals
selDs           = selBasic + " && (gen_sig == 0 || gen_sig == 1)"                                  # select Ds signals
selDsStar       = selBasic + " && (gen_sig == 5 || gen_sig == 6)"                                  # select Ds star signals

## BASELINE

baseline = ' & '.join([
f'(dsMu_m < {bsMass_})',
'(k1_charge*k2_charge <0)',
'(mu_charge*pi_charge < 0)',
'(mu_pt > 8)',
'(k1_pt > 1)',
'(k2_pt > 1)',
'(pi_pt > 1)',
'(lxy_ds < 1)',
'(mu_id_medium == 1)',
'(mu_rel_iso_03 < 0.3)',
#'tv_prob > 0.1',
#'((cosPiK1 < -0.3) || (cosPiK1 > 0.3))',
'(fv_prob > 0.1)'
])

baselineHb =        baseline + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 5) && (gen_sig != 6) & (gen_match_success ==1 ) "
baselineDsMu =      baseline + " && (gen_sig == 0)" 
baselineDsTau =     baseline + " && (gen_sig == 1)"
baselineDsStarMu =  baseline + " && (gen_sig == 5)" 
baselineDsStarTau = baseline + " && (gen_sig == 6)"

#########################################
## CREATE RDF FROM TREE                ##
#########################################

def getRdf(dateTime, debug = False):

  #access the flat ntuples
  files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test

  if debug:
    print("picking only one file for debugging ...")
    fileList = os.listdir(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/")
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/" +  fileList[0] # pick just the first one 

  #chain them
  chain = ROOT.TChain("tree")
  chain.Add(files)

  #create rdf from tree
  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)


# Create rdf from tree
print(f" ===> Start creating RDataFrames")
chainSig, rdfSig   = getRdf(sig)#, debug = True)
chainHb,  rdfHb    = getRdf(hb)#, debug = True)
chainData,rdfData  = getRdf(data)#, debug = True)

#########################################
## COLOR SCHEME                        ##
#########################################


def getColorAndLabel(var):

  "Colors for simple plot of signgle histograms, default" 

  if   "gen"      in var: color,label = ROOT.kMagenta,    "Gen"
  elif "coll"     in var: color,label = ROOT.kBlack,      "Coll."
  elif "lhcb_alt" in var: color,label = ROOT.kBlue,       "LHCb xyz"
  elif "lhcb"     in var: color,label = ROOT.kViolet,     "LHCb z"
  elif "reco_1"   in var: color,label = ROOT.kGreen,      "Math. 1"
  elif "reco_2"   in var: color,label = ROOT.kOrange,     "Math. 2"
  elif "weighted" in var: color,label = ROOT.kRed,        "Weighted Math."

  else:                   color,label = ROOT.kBlack,      "Reco"

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
 
  elif "stacked" in key:
    colors = {}
   
    labels.append(r"D_{s} #mu"); labels.append(r"D_{s} #tau"); labels.append(r"D*_{s} #mu");   labels.append(r"D*_{s} #tau"); labels.append(r"H_{b}"); labels.append(r"Data");
    colors = {"dsMu":   ROOT.kAzure,
              "dsTau":  ROOT.kPink,
              "dsStarMu":   ROOT.kMagenta - 4,
              "dsStarTau":  ROOT.kGreen + 3,
              "hb":         ROOT.kRed,
              "combL":      ROOT.kGray,
              "combR":      ROOT.kGray,
              "data":       ROOT.kBlack}

  if "hb" in key:
    hbColor,hbLabel = getHbColorAndLabel()
    labels.append(hbLabel)
    colors.append(hbColor)

  return (colors,labels)

#########################################
## CREATE LEGEND                       ##
#########################################

def getLegend(key, x0 = 0.5, y0 = 0.8, x1 = 0.9, y1 = 0.9, nColumns = 1):

  "Default setting is x0,y0,x1,y1 as given with 1 column"

  if ('all' in key or 'method' in key):
    # 4 signals and no gini, stretch legend
    x0 = 0.3
    y0 = 0.75
    

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


def getReco1Scale(name, selection, start, stop, chain = chainSig):

  ntot   = chain.GetEntries(selection + f"&&      ({name}_reco_1 > {start}) && ({name}_reco_1 < {stop}) && (gen_{name} > {start}) && (gen_{name} < {stop}) && ({name}_reco_1 == {name}_reco_1)")
  nReco1 = chain.GetEntries(selection + f"&& (abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && ({name}_reco_1 == {name}_reco_1)")
    
  return nReco1 / ntot 

def getWeightedReco(histos, name, selection, norm = True):

  #get bin edges
  nBins = histos[name + "_reco_1"].GetNbinsX()
  start = histos[name + "_reco_1"].GetXaxis().GetBinLowEdge(1) #first bin ( 0 is underflow!)
  stop  = histos[name + "_reco_1"].GetXaxis().GetBinUpEdge(nBins) #last bin ( 0 is underflow!)
  print(start)
  print(stop)
  #get reco scale
  reco1Scale = getReco1Scale(name, selection, start, stop)

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

def createHistos(selection,rdf, linewidth = 2):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"

  histos = {}
   
  for var,model in models.items(): 

      histos[var] = rdf.Filter(selection).Histo1D(model[0], var)
      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetMaximum())

      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(linewidth)
  print(color)
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
  legend = getLegend("simple", nColumns = 2)
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
  print(key)
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend(key, nColumns = 2)
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
  legend = getLegend("method",nColumns = 3)
  canvas.cd()

  yMax = getYMax([histos[var] for var in list(histos.keys())],norm)

  #recreate the weighted combination
  recoAvailable = name + "_reco_1" in models.keys() #check that recos are available 
    
  for i, var in enumerate(list(histos.keys())):

    hist = histos[var]

    if norm: 
      hist.Scale(1/hist.Integral())
      hist.GetYaxis().SetTitle("a.u.")

    hist.SetMaximum(1.3*yMax)    
    _,label = getColorAndLabel(var)

    # draw 
    if (i == 0): hist.Draw("HIST")
    else:        hist.Draw("HIST SAME")

    legend.AddEntry(hist.GetValue(),label,"l")
    
  if recoAvailable: 
    weightedReco = getWeightedReco(histos,name,selection,norm)
    
    weightedReco.Draw("HIST SAME") 
    weightedReco.SetMaximum(1.3*yMax)
    col,label = getColorAndLabel("weighted_reco")
    weightedReco.SetLineColor(col)

    legend.AddEntry(weightedReco,label,"l")

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
#histosBasicHb       = createHistos(selBasicHb,  rdfHb)

print(f" ===> Basic, Mu, Tau, Ds, DsStar DONE")

#histosDsMu      = createHistos(selDsMu,      rdfSig)
#histosDsTau     = createHistos(selDsTau,     rdfSig)
#histosDsStarMu  = createHistos(selDsStarMu,  rdfSig)
#histosDsStarTau = createHistos(selDsStarTau, rdfSig)
print(f" ===> Single Signals DONE")

# Simple Plots 
#callSimpleDraw(histosBasic,histosBasicHb, selBasic,selBasicHb)
"""
# Plot Mu vs Tau
callSignalDistinction( [histosMu, histosTau, histosBasicHb],    [selBasic,selBasicHb], "mutau")
# Plot Ds vs Ds star
callSignalDistinction( [histosDs, histosDsStar, histosBasicHb], [selBasic,selBasicHb], "dsdsstar")
# Plot all signals
callSignalDistinction( [histosDsMu, histosDsTau, histosDsStarMu, histosDsStarTau, histosBasicHb], [selBasic,selBasicHb], "allsignals")

#Method comparison
prefix = ["q2","m2_miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt", "e_star"]
for name in prefix:
  print(f" ===> Compare Bs reco methods for {name}")
  # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...
  histos = {var: histosDsMu[var] for var in models.keys() if name in var}
  methodDistinction(histos,name,selDsMu)

#########################################
## Prepare stacked plots (do SB fit)    #
#########################################
"""

def getHbScale(selHb, selMu):

  # get plotting limits of phiPi mass

  bins  = models["phiPi_m"][0].fNbinsX
  start = models["phiPi_m"][0].fXLow
  stop  = models["phiPi_m"][0].fXUp

  pureHb = ROOT.TH1D("pureHb","pureHb",bins,start,stop)
  pureMu = ROOT.TH1D("pureMu","pureMu",bins,start,stop)
  muNew = ROOT.TH1D("muNew","muNew",bins,start,stop)

  chainHb.Project("pureHb", "phiPi_m", selHb   + f" & (phiPi_m > {start} ) & (phiPi_m < {stop} )" )
  chainHb.Project("pureMu", "phiPi_m", selMu + f" & (phiPi_m > {start} ) & (phiPi_m < {stop} )")
  chainSig.Project("muNew",    "phiPi_m", selMu + f" & (phiPi_m > {start} ) & (phiPi_m < {stop} )" )

  scale_hb = muNew.Integral() * pureHb.Integral() / pureMu.Integral()
  print( "#Hb events for MC sample is:", scale_hb)
  return scale_hb



sigma, h          = getSigma(rdfSig, "phiPi_m", baseline + "& gen_sig == 0")
#A, B, C, S        = getABCS( data, rdfData, selection , "phiPi_m", sigma, h, bins = 15)
A = 24158.821588264094; B = 226543.39752377832; C = 21219.183790935822; S = 19829.186060030886; #use for debugging (faster) 

tresh = n * sigma
mlow = dsMass_ - tresh
mlow2 = mlow - 1*sigma
mhigh = dsMass_ + tresh
mhigh2 = mhigh + 1*sigma

signalRegion = f"& ({mlow} < phiPi_m) & (phiPi_m < {mhigh})"
leftSB       = f"& ({mlow2} < phiPi_m) & (phiPi_m < {mlow})"
rightSB      = f"& ({mhigh} < phiPi_m) & (phiPi_m < {mhigh2})"

# get fakemass histo
fakemass = np.genfromtxt('mass_bincontent.csv', delimiter=',')

#create histos returns a dictionary !:)

#for all variables expect mass plot only signal region (indicated with 'S')
baselineSDsMu      = createHistos(baselineDsMu      + signalRegion, rdfSig)
baselineSDsTau     = createHistos(baselineDsTau     + signalRegion, rdfSig)
baselineSDsStarMu  = createHistos(baselineDsStarMu  + signalRegion, rdfSig)
baselineSDsStarTau = createHistos(baselineDsStarTau + signalRegion, rdfSig)

baselineSHb        = createHistos(baselineHb        + signalRegion, rdfHb)
baselineDataL      = createHistos(baseline          + leftSB, rdfData)
baselineDataR      = createHistos(baseline          + rightSB, rdfData)
baselineSData       = createHistos(baseline         + signalRegion, rdfData)

# for the ds mass plot we also want to plot the sidebands! (indicated with 'C' for complete)
baselineCDsMu      = createHistos(baselineDsMu      , rdfSig)
baselineCDsTau     = createHistos(baselineDsTau     , rdfSig)
baselineCDsStarMu  = createHistos(baselineDsStarMu  , rdfSig)
baselineCDsStarTau = createHistos(baselineDsStarTau , rdfSig)

baselineCHb        = createHistos(baselineHb        , rdfHb)
baselineCData       = createHistos(baseline         , rdfData)

hb_scale = getHbScale(baselineHb, baselineDsMu)

def stackedPlot(histos, var, hb_scale):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    histos[key].GetValue().SetFillColor(color[key])
    histos[key].GetValue().SetLineColor(color[key])
    if key == "data":
       histos[key].GetValue().SetMarkerStyle(8) 

  # take as model
  bins  = histos["dsMu"].GetNbinsX()
  start = histos["dsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  stop  = histos["dsMu"].GetXaxis().GetBinUpEdge(bins) # last bin


  # stacked histo
  hs = ROOT.THStack(var,"")
  # error histo
  hErr = ROOT.TH1D(var,"",bins,start,stop)
  hErr.SetName("err")
  hErr.GetXaxis().SetTitleOffset(1.3)
  hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())  
  hErr.SetLineWidth(0)
  # mass histo for fakemass
  hComb = hErr.Clone()
  hComb.SetName("mass")
  hComb.SetFillColor(ROOT.kGray)
  hComb.SetLineWidth(0)
  # special error fillstyle 
  hErr.SetFillColor(ROOT.kGray+1)
  hErr.SetFillStyle(3344)

  #scale hb to other mc

  print("Scale of Hb: ", hb_scale / histos["hb"].Integral())
  histos["hb"].Scale(hb_scale / histos["hb"].Integral())

  nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()
  print("nsig = ", nSig)

  ###################################
  ## Combinatorial treatment       ##
  ###################################

  if var != "phiPi_m":
    print("i am here")
    # use sideband method
    # left one
    hComb = histos["combL"].Clone()
    hComb.SetName("combS") #signal region
    hComb.Scale(B/(2*A))

    print("i am here")
    # right one
    hCombR = histos["combR"].Clone()
    hCombR.Scale(B/(2*C))
    hComb.Add(hCombR)

    print("i am here")
    nComb = hComb.Integral()

    print("i am here")
    # add comb to stack and err 
    hErr.Add(hComb)
    hs.Add(hComb)
 
  else:
    for i in range(bins):
      hComb.SetBinContent(i+1,fakemass[i])
      hComb.SetBinError(i+1,np.sqrt(fakemass[i]))
      print( "data now:", histos["data"].GetBinContent(i+1))
      print("fakemass now:", hComb.GetBinContent(i+1))
    
    hComb.Scale(0.999)

    nComb = hComb.Integral()

    # add comb to stack and err 
    hErr.Add(hComb)
    hs.Add(hComb)

  ###################################
  ## Scale MC to data              ##
  ###################################

  for key in histos.keys():
    if ('comb' not in key) and ('data' not in key):

      histos[key].Scale((histos["data"].Integral() - nComb) / nSig)
      if var == 'phiPi_m': histos[key].Scale(0.999)
      hErr.Add(histos[key].GetValue())
      hs.Add(histos[key].GetValue())
 
  ###################################
  ## Create Pads                   ## 
  ###################################

  c1 = ROOT.TCanvas(var, '', 700, 700)
  c1.Draw()
  c1.cd()
  main_pad = ROOT.TPad('main_pad', '', 0., 0.25, 1. , 1.  )
  main_pad.Draw()
  c1.cd()
  ratio_pad = ROOT.TPad('ratio_pad', '', 0., 0., 1., 0.25)
  ratio_pad.Draw()
  main_pad.SetTicks(True)
  main_pad.SetBottomMargin(0.)
  main_pad.SetLeftMargin(.16)
  ratio_pad.SetTopMargin(0.)
  ratio_pad.SetLeftMargin(.16)
  ratio_pad.SetGridy()
  ratio_pad.SetBottomMargin(0.45)

  main_pad.cd()

  ###################################
  ## legends                       ## 
  ###################################
 
  #legend
  leg = ROOT.TLegend(.2,.72,.88,0.88)
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(42)
  leg.SetTextSize(0.035)
  leg.SetNColumns(4)
  leg.AddEntry(histos["dsStarMu"].GetValue()   ,'B_{S}#rightarrow D*_{s}#mu#nu' ,'F' )
  leg.AddEntry(histos["dsMu"].GetValue()       ,'B_{S}#rightarrow D_{s}#mu#nu'  ,'F' )
  leg.AddEntry(histos["hb"].GetValue()         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  leg.AddEntry(histos["data"].GetValue()       ,'Data','LEP')
  leg.AddEntry(histos["dsStarTau"].GetValue()  ,'B_{S}#rightarrow D*_{s}#tau#nu','F' )
  leg.AddEntry(histos["dsTau"].GetValue()      ,'B_{S}#rightarrow D_{s}#tau#nu' ,'F' )
  leg.AddEntry(hComb                           ,'Comb. Bkg.'  ,'F' )
  leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
  #plot mainpad
  histos["data"].GetYaxis().SetRangeUser(1e-3, histos["data"].GetMaximum()*1.6)
  histos["data"].GetYaxis().SetTitle('Events')
  
  #draw data with uncertainty
  histos["data"].Draw("EP")
  #draw stacked histo
  hs.Draw("HIST SAME")
  #draw data again (used for plotting uses)
  histos["data"].Draw("EP SAME")
  #draw stat uncertainty of MC, (can not be directly done with a stacked histo)
  hErr.Draw("E2 SAME")
  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = 'X fb^{-1}')
  
  #plot ratiopad
  ratio_pad.cd()
  ratio = histos["data"].Clone()
  ratio.SetName(ratio.GetName()+'_ratio')
  ratio.Divide(hErr)
  ratio_stats = hErr.Clone()
  ratio_stats.SetName(ratio.GetName()+'_ratiostats')
  ratio_stats.Divide(hErr)
  ratio_stats.SetMaximum(1.999) # avoid displaying 2, that overlaps with 0 in the main_pad
  ratio_stats.SetMinimum(0.0001) # and this is for symmetry
  ratio_stats.GetYaxis().SetTitle('Data / MC')
  ratio_stats.GetYaxis().SetTitleOffset(0.5)
  ratio_stats.GetYaxis().SetNdivisions(405)
  ratio_stats.GetXaxis().SetLabelSize(3.* ratio.GetXaxis().GetLabelSize())
  ratio_stats.GetYaxis().SetLabelSize(3.* ratio.GetYaxis().GetLabelSize())
  ratio_stats.GetXaxis().SetTitleSize(3.* ratio.GetXaxis().GetTitleSize())
  
         
  #divide signals by total number of measurements
  norm_stack = ROOT.THStack('norm_stack', '')
 
  for key in histos.keys():
    if ('comb' not in key) and ('data' not in key):
      h = histos[key].Clone()
      h.Divide(hErr)
      norm_stack.Add(h)

  hDummy = hComb.Clone() #do it also for the comb
  hDummy.Divide(hErr)
  norm_stack.Add(hDummy)
  
  norm_stack.Draw('hist same')
  
  #subplot
  line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
  line.SetLineColor(ROOT.kBlack)
  line.SetLineWidth(1)
  ratio_stats.GetYaxis().CenterTitle()
  
  ratio_stats.Draw('EP')
  norm_stack.Draw('hist same')
  ratio_stats.Draw('EP same')
  line.Draw('same')
  ratio.Draw('EP same')
  ROOT.gPad.RedrawAxis()
  
  #saving
  c1.Modified()
  c1.Update()
  c1.SaveAs(var + ".pdf")
  print(f"===> Produced plot: {var}.pdf")


for var in models.keys():
  print(f"===> Producing stacked plot for variable: {var}") 
  if var != "phiPi_m":

    histos = {"dsMu":      baselineSDsMu[var], 
             "dsTau":     baselineSDsTau[var],
             "dsStarMu":  baselineSDsStarMu[var], 
             "dsStarTau": baselineSDsStarTau[var], 
             "hb"       : baselineSHb[var], 
             "combL"    : baselineDataL[var], 
             "combR"    : baselineDataR[var], 
             "data"     : baselineSData[var]}
    stackedPlot(histos, var, hb_scale)

  else:

    histos = {"dsMu":      baselineCDsMu[var], 
             "dsTau":     baselineCDsTau[var],
             "dsStarMu":  baselineCDsStarMu[var], 
             "dsStarTau": baselineCDsStarTau[var], 
             "hb"       : baselineCHb[var], 
             "data"     : baselineCData[var]}
    stackedPlot(histos, var, hb_scale)

