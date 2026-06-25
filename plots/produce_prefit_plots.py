import glob
import ROOT
import argparse
import os
import re
import sys
import yaml
import uproot
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
ROOT.ROOT.EnableImplicitMT() 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio, getSignflipRatioTest, fitAnotherVar
from helper import * 
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels
from blinding import *

import numpy as np
from cms_style import CMS_lumi

#ROOT.ROOT.EnableImplicitMT(8)

def boolean_string(s):
    if s not in {"false", "true"}:
        raise ValueError("Error: Not a valid boolean string, please use 'true' or 'false' (all lowercase!)")
    return s == "true"

# parsing
parser = argparse.ArgumentParser()
parser.add_argument("--datetime",   help = "") 
parser.add_argument("--label"   ,   help = "") 
args = parser.parse_args()
dt = args.datetime

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

#where to save plots and datacards
toSave_plots        = f"/work/pahwagne/RDsTools/plots/cmsplots_binned/{dt}/"
toSave_cards        = f"/work/pahwagne/RDsTools/fit/datacards_binned/{dt}"
shapes_folder       = f"/work/pahwagne/RDsTools/fit/shapes_binned/{dt}/"
shapes_folder_blind = f"/work/pahwagne/RDsTools/fit/shapes_binned/{dt}/blind/"

if not os.path.exists(toSave_plots): 
  os.makedirs(toSave_plots)
  os.makedirs(toSave_plots + "/log")  

if not os.path.exists(shapes_folder): 
  os.makedirs(shapes_folder)

if not os.path.exists(shapes_folder_blind): 
  os.makedirs(shapes_folder_blind)

# load systematics
print(f"====> Adding the following systematics: {systematics_scalar} and {systematics_vector}") #defined in helper.py
sys_scalar = systematics_scalar
sys_vector = systematics_vector


massYields = {}
modifyHammer = False

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
    colors = {"DsMu":            ROOT.TColor.GetColor("#92dadd"), #ROOT.kBlue - 2,
              "DsMu_woHammer":   ROOT.TColor.GetColor("#92dadd"), #ROOT.kBlue - 2,
              "Mu_in_Hb":        ROOT.TColor.GetColor("#92dadd"), #ROOT.kBlue - 2,
              "DsTau":           ROOT.TColor.GetColor("#ffa90e"), #ROOT.kGreen,
              "DsStarMu":        ROOT.TColor.GetColor("#3f90da"), #ROOT.kCyan,
              "DsStarTau":       ROOT.TColor.GetColor("#e76300"), #ROOT.kOrange,

              "Hb":              ROOT.TColor.GetColor("#bd1f01"), #ROOT.kRed,
              "Hb_fd":           ROOT.TColor.GetColor("#bd1f01"), #ROOT.kRed -7,
              "Hb_dc":           ROOT.TColor.GetColor("#bd1f01"), #ROOT.kRed -2,

              "Hb_others":       ROOT.TColor.GetColor("#bd1f01"), #ROOT.kRed -5,

              "Hb_bs":           ROOT.TColor.GetColor("#bd1f01"), #ROOT.kRed +2,
              "Hb_bs_fd":        ROOT.TColor.GetColor("#832db6"), #ROOT.kRed +2,
              "Hb_bs_dc":        ROOT.TColor.GetColor("#a96b59"), #ROOT.kRed -7,

              "Hb_b0":           ROOT.TColor.GetColor("#bd1f01"), #ROOT.kMagenta,
              "Hb_b0_fd":        ROOT.TColor.GetColor("#b9ac70"), #ROOT.kMagenta,
              "Hb_b0_dc":        ROOT.TColor.GetColor("#94a4a2"), #ROOT.kMagenta -7,

              "Hb_bpm":          ROOT.TColor.GetColor("#bd1f01"), #ROOT.kOrange + 7,
              "Hb_bpm_fd":       ROOT.TColor.GetColor("#e42536"), #ROOT.kOrange + 7,
              "Hb_bpm_dc":       ROOT.TColor.GetColor("#964a8b"), #ROOT.kOrange + 5,

              "Hb_lambdab":      ROOT.TColor.GetColor("#bd1f01"), #ROOT.kViolet,
              "Hb_lambdab_fd":   ROOT.TColor.GetColor("#7a21dd"), #ROOT.kViolet,
              "Hb_lambdab_dc":   ROOT.TColor.GetColor("#9c9ca1"), #ROOT.kViolet-7,

              "bs":              ROOT.TColor.GetColor("#3f90da"), #ROOT.kRed - 9,
              "b0":              ROOT.TColor.GetColor("#3f90da"), #ROOT.kRed - 6,
              "bplus":           ROOT.TColor.GetColor("#3f90da"), #ROOT.kRed - 5,
              "combL":           ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "combR":           ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data_sf":         ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data_sf_kk":      ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data_sf_pimu":    ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data_sf_both":    ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data_sf_one":     ROOT.TColor.GetColor("#717581"), #ROOT.kGray+1,
              "Data":            ROOT.kBlack}

  if "Hb" in key:
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
    if norm: 
      try: hist.Scale(1/hist.Integral())
      except: pass 
    maxis.append(hist.GetMaximum())
 
  return max(maxis)

#########################################
## TOOLS FOR WEIGHTED RECO METHOD      ##
#########################################



def getWeightedReco(histos, name, selection, norm = True):

  #get bin edges
  nBins = histos[name + "_reco_1"].GetNbinsX()
  start = histos[name + "_reco_1"].GetXaxis().GetBinLowEdge(1) #first bin ( 0 is underflow!)
  stop  = histos[name + "_reco_1"].GetXaxis().GetBinUpEdge(nBins) #last bin ( 0 is underflow!)
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
## Lambda functions for hammer weight  ##
#########################################

ROOT.gROOT.ProcessLine(r'''float get_hammer_weight(float central, int sig){

  float average_w;
  float w = central;

  double dsStarTau_w = 3.304854088595039;
  double dsStarMu_w  = 3.198165968764498;

  if (sig == 0){ 
    w = 1.0;
    average_w = 1.0; 
  }
  else if (sig == 1){ 
    w = 1.0;
    average_w = 1.0; 
  }
  else if (sig == 10){ 
    average_w = dsStarMu_w;
  }
  else if (sig == 11){ 
    average_w = dsStarTau_w;
  }

  return w / average_w;

  }''')


###########################################
# Prepare signflip fit                    #
###########################################

def prepareSignFlip (hb, hb_scale, signals):

  #scale hb to other MC
  if hb_scale == 0:
    hb.Scale(0.0 )
  else:
    hb.Scale(hb_scale / hb.Integral())

  #add all signals and hb together
  hRest    = hb 
  for sig in signals:
    hRest.Add( sig ) 

  return hRest
 
#def addSystematics(hist_dict, var, selec_DsTau, selec_DsMu, selec_DsStarTau, selec_DsStarMu ):
#
#   # if there are systematics, include them in the scaling
#   for s in sys_scalar:
#     for direc in ["Up", "Down"]:
#       hist_dict["dsTau_"     + s + direc] = selec_DsTau    [var + "_" + s + direc]
#       hist_dict["dsMu_"      + s + direc] = selec_DsMu     [var + "_" + s + direc]
#
#   for s in sys_vector:
#     for direc in ["Up", "Down"]:
#       hist_dict["dsStarTau_" + s + direc] = selec_DsStarTau[var + "_" + s + direc]
#       hist_dict["dsStarMu_"  + s + direc] = selec_DsStarMu [var + "_" + s + direc]
#
#   return hist_dict


#########################################
## CREATE DEFAULT HISTOS               ##
#########################################



def plotSignFlipShapes(kk,pimu,both,var):
    c = ROOT.TCanvas("test", "test", 800,800)

    print("Flip shape for variable:", var)
    h1 = kk[var].Clone() 
    h1.Scale(1/h1.Integral())
    h1.SetLineColor(ROOT.kRed)

    h2 = pimu[var].Clone()
    h2.Scale(1/h2.Integral())
    h2.SetLineColor(ROOT.kBlue)

    h3 = both[var].Clone()
    h3.Scale(1/h3.Integral())
    h3.SetLineColor(ROOT.kGreen)

    lege = getLegend("easy")
    lege.AddEntry(h1, "KK wrong", "L")
    lege.AddEntry(h2, "#pi #mu wrong", "L")
    lege.AddEntry(h3, "KK and #pi#mu wrong", "L")
    h1.SetMaximum(0.20)
    h1.Draw("NORM HIST")
    h2.Draw("NORM HIST SAME")
    h3.Draw("NORM HIST SAME")
    lege.Draw("SAME")
    c.SaveAs(f"shapes_wrong_sign_{var}.pdf") 



def stackedPlot(histos2, var, scale_hb, scale_pimu = None, scale_rest = None, logx = False, logy = False, region = None, blind = False):

  histos = {}
  for key, hist in histos2.items():

    #print(key) 
    # we skip the blind keys and rename the DsTau_blind to DsTau for the corresponding cases
    if "blind" in key: continue

    if   (key == "DsTau"     and blind == True) : histos[key] = histos2["DsTau_blind"].Clone()
    elif (key == "DsStarTau" and blind == True) : histos[key] = histos2["DsStarTau_blind"].Clone()
    else                                        : histos[key] = hist.Clone()
    #print(histos) 

  #print("==========> number of hb events right calling stackedPlot", histos["Hb"].Integral())
  #print("==========> number of bs events right calling stackedPlot", histos["Hb_bs"].Integral())
  #print("==========> number of b0 events right calling stackedPlot", histos["Hb_b0"].Integral())
  #print("==========> number of bpm events right calling stackedPlot", histos["Hb_bpm"].Integral())
  #print("==========> number of bpm events right calling stackedPlot", histos["Hb_lambdab"].Integral())
  #print("==========> number of other events right calling stackedPlot", histos["Hb_others"].Integral())

  color, labels = getColorAndLabelSignalDistinction("stacked")

  #supress data plotting in SR. region is a string of type: ch0
  dataBlind = False
  if ((int(region[2:]) > 5) and (int(region[2:]) < 30)): 
    dataBlind = True 
    #color["DsTau"] = color["DsMu"]
    #color["DsStarTau"] = color["DsStarMu"]

  # if we are anyway using the blinded histos (with a blind scale, no need do hide data)
  if (blind == True): dataBlind = False

  for i,key in enumerate(histos.keys()):

    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillColor(color[key])
    histos[key].SetLineColor(color[key])
    if key == "Data":
       histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["DsMu"].GetNbinsX()

  # stacked histo
  hs = ROOT.THStack(var,"")

  # error histo
  hErr = histos["DsMu"].Clone()  
  hErr.Reset()
  hErr.SetName("err")
  hErr.GetXaxis().SetTitleOffset(1.3)
  hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())  
  hErr.SetLineWidth(0)
  hErr.SetFillColor(ROOT.kBlack)
  hErr.SetFillStyle(3144)
  hErr.SetMarkerStyle(0)

  ###################################
  ## scale Hb to signals           ##
  ###################################

  for key in histos.keys():
    if ("Hb" in key) : 
      #print(f"===> scale {key} to signals")
      histos[key].Scale(scale_hb)

  #nSig = histos["DsMu"].Integral() + histos["DsTau"].Integral() + histos["DsStarMu"].Integral() + histos["DsStarTau"].Integral() + histos["Hb"].Integral()
  #print("==========> number of tau* signals before scaling", histos["DsStarTau"].Integral())
  #print("==========> number of mu* signals before scaling", histos["DsStarMu"].Integral())
  #print("==========> number of hb events after first sclaing", histos["Hb"].Integral())
  #print("==========> number of bs events after first sclaing", histos["Hb_bs"].Integral())
  #print("==========> number of b0 events after first sclaing", histos["Hb_b0"].Integral())
  #print("==========> number of bpm events after first sclaing", histos["Hb_bpm"].Integral())
  #print("==========> number of bpm events right calling stackedPlot", histos["Hb_lambdab"].Integral())
  #print("==========> number of other event right calling stackedPlot", histos["Hb_others"].Integral())


  ###################################
  ## Scale Comb to data            ##
  ###################################

  # normalize
 
  hComb = histos["Data_sf_pimu"].Clone()
  hComb.Scale(scale_pimu)

  #also scale the Up/Down
  for key in histos.keys():
    if "Data_sf_pimu" in key:
      histos[key].Scale(scale_pimu)

  hComb.SetLineColor(ROOT.kGray + 1)
  hComb.SetFillColor(ROOT.kGray + 1)

  nComb = hComb.Integral()
  histos["comb"] = hComb # append to create datacard later

  ###################################
  ## Scale MC to data              ##
  ###################################

  for key in histos.keys():

    #normalize and scale the signals and hb to data
    if ('comb' not in key) and ('Data' not in key):
      #print(f"===> sclaing {key} to data")
      histos[key].Scale(scale_rest) # NEW

  ###################################
  ## Scale everythin to data       ##
  ###################################

  #ntot = histos["DsMu"].Integral() + histos["DsTau"].Integral() + histos["DsStarMu"].Integral() + histos["DsStarTau"].Integral() + histos["Hb"].Integral() + histos["comb"].Integral() 
  #for key in histos.keys():
  #  if 'Data' not in key:
  #    histos[key].Scale(histos["Data"].Integral() / ntot )

  # add comb to stack and err 

  ##################################
  ## Build stacked (err) Histos   ##
  ##################################

  #explicit order vor visibility
  hs.Add(hComb)
  hs.Add(histos["DsTau"     ])
  hs.Add(histos["DsMu"      ])
  hs.Add(histos["DsStarTau" ])
  hs.Add(histos["DsStarMu"  ])
  #hs.Add(histos["Hb"     ])
  hs.Add(histos["Hb_bs_fd"     ])
  hs.Add(histos["Hb_bs_dc"     ])

  hs.Add(histos["Hb_b0_fd"     ])
  hs.Add(histos["Hb_b0_dc"     ])

  hs.Add(histos["Hb_bpm_fd"    ])
  hs.Add(histos["Hb_bpm_dc"    ])

  hs.Add(histos["Hb_lambdab_fd"])
  hs.Add(histos["Hb_lambdab_dc"])

  hs.Add(histos["Hb_others" ])

  #(here order doesnt matter)
  hErr.Add(hComb)
  hErr.Add(histos["DsTau"     ])
  hErr.Add(histos["DsStarTau" ])
  hErr.Add(histos["DsMu"      ])
  hErr.Add(histos["DsStarMu"  ])
  #hErr.Add(histos["Hb"     ])
  hErr.Add(histos["Hb_bs_fd"     ])
  hErr.Add(histos["Hb_bs_dc"     ])

  hErr.Add(histos["Hb_b0_fd"     ])
  hErr.Add(histos["Hb_b0_dc"     ])

  hErr.Add(histos["Hb_bpm_fd"    ])
  hErr.Add(histos["Hb_bpm_dc"    ])

  hErr.Add(histos["Hb_lambdab_fd"])
  hErr.Add(histos["Hb_lambdab_dc"])

  hErr.Add(histos["Hb_others" ])


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
  leg = ROOT.TLegend(.2,.60,.88,0.88)
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(42)
  leg.SetTextSize(0.035)
  leg.SetNColumns(2)
  leg.AddEntry(histos["DsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'F' )
  leg.AddEntry(histos["DsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'F' )
  leg.AddEntry(histos["DsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','F' )
  leg.AddEntry(histos["DsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'F' )

  #leg.AddEntry(histos["Hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  #leg.AddEntry(histos["Hb_fd"]       ,'B^{#pm 0}_{(s)}, #Lambda_{b} #rightarrow feed-down'   ,'F' )
  #leg.AddEntry(histos["Hb_dc"]       ,'B^{#pm 0}_{(s)}, #Lambda_{b} #rightarrow double-charm','F' )
  leg.AddEntry(histos["Hb_bs_fd"]       ,'B_{s} #rightarrow feed-down','F' )
  leg.AddEntry(histos["Hb_bs_dc"]       ,'B_{s} #rightarrow double-charm','F' )

  leg.AddEntry(histos["Hb_bpm_fd"]      ,'B^{#pm} #rightarrow feed-down','F' )
  leg.AddEntry(histos["Hb_bpm_dc"]      ,'B^{#pm} #rightarrow double-charm','F' )

  leg.AddEntry(histos["Hb_b0_fd"]       ,'B^{0} #rightarrow feed-down','F' )
  leg.AddEntry(histos["Hb_b0_dc"]       ,'B^{0} #rightarrow double-charm','F' )

  leg.AddEntry(histos["Hb_lambdab_fd"]  ,'#Lambda_{b} #rightarrow feed-down','F' )
  leg.AddEntry(histos["Hb_lambdab_dc"]  ,'#Lambda_{b} #rightarrow double-charm','F' )
 
  leg.AddEntry(histos["Hb_others"]   ,'other b #rightarrow D_{s} + #mu','F' )
  leg.AddEntry(histos["Data"]        ,'Data','LEP')
  #leg.AddEntry(histos["comb"]                           ,'Comb. + Fakes'  ,'F' )
  leg.AddEntry(hComb                            ,'Comb. + Fakes'  ,'F' )
  #leg.AddEntry(hComb_kk                             ,'kk'  ,'F' )
  leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
  #plot mainpad
  if logy:
    ROOT.gPad.SetLogy()
    histos["Data"].GetYaxis().SetTitle('Events')
    histos["Data"].GetYaxis().SetRangeUser(1e-3, hs.GetMaximum()*1000)
    histos["Data"].SetMinimum(100.001) # avoid idsplaying tons of comb bkg

  elif logx:
    ROOT.gPad.SetLogx()

    yAxisTitle = "events"
    histos["Data"].GetYaxis().SetTitle(yAxisTitle)
    histos["Data"].GetYaxis().SetRangeUser(1e-3, histos["Data"].GetBinContent(histos["Data"].GetMaximumBin())*1.5)

  else:

    yAxisTitle = "events"
    histos["Data"].GetYaxis().SetTitle(yAxisTitle)
    histos["Data"].GetYaxis().SetRangeUser(1e-3, histos["Data"].GetBinContent(histos["Data"].GetMaximumBin())*1.8)

  #if we are in the signal region and do not blind the MC, then we blind the data!

  if dataBlind and not blind:
    #only draws axis, but not the Data! :D
    histos["Data"].Draw("AXIS")

  else:
    #draw Data with uncertainty
    histos["Data"].Draw("EP")
  #histos["Data"].Draw("EP")

  #draw stacked histo
  hs.Draw("HIST SAME")

  #draw Data again (used for plotting uses)
  if dataBlind and not blind:
    histos["Data"].Draw("AXIS SAME")
  else:
    histos["Data"].Draw("EP SAME")
  #histos["Data"].Draw("EP SAME")


  #draw stat uncertainty of MC, (can not be directly done with a stacked histo)
  hErr.Draw("E2 SAME")
  #hErr.Draw("E2 SAME")

  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  
  #if trigger == "mu7": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Private Work', lumi_13TeV = '6.4 fb^{-1}')
  #if trigger == "mu9": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Private Work', lumi_13TeV = '20.7 fb^{-1}')
  #CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Private Work', lumi_13TeV = '2018D, Part 1')
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Private Work', lumi_13TeV = f'{args.label}')
 
 
  #plot ratiopad
  ratio_pad.cd()

  #ROOT.gPad.SetLogy()
  ratio = histos["Data"].Clone()
  ratio.SetName(ratio.GetName()+'_ratio')
  ratio.Divide(hErr)
  ratio_stats = hErr.Clone()
  ratio_stats.SetName(ratio.GetName()+'_ratiostats')
  ratio_stats.Divide(hErr)
  ratio_stats.SetMaximum(1.299) # avoid displaying 2, that overlaps with 0 in the main_pad
  ratio_stats.SetMinimum(0.7) # and this is for symmetry
  ratio_stats.GetYaxis().SetTitle('Data / MC')
  ratio_stats.GetXaxis().SetTitle(ratio.GetXaxis().GetTitle())

  ratio_stats.GetYaxis().SetTitleOffset(0.5)
  ratio_stats.GetYaxis().SetNdivisions(405)
  ratio_stats.GetYaxis().SetTitleOffset(0.5)

  ratio_stats.GetXaxis().SetLabelSize(3.* ratio.GetXaxis().GetLabelSize())
  ratio_stats.GetYaxis().SetLabelSize(3.* ratio.GetYaxis().GetLabelSize())
  ratio_stats.GetXaxis().SetTitleSize(3.* ratio.GetXaxis().GetTitleSize())
  ratio_stats.GetYaxis().SetTitleSize(3.* ratio.GetYaxis().GetTitleSize())
  
         
  #divide signals by total number of measurements
  norm_stack = ROOT.THStack('norm_stack', '')
 
  #explicit order for visibility
  for key in [
    "comb",
    "DsStarMu",
    "DsStarTau",
    "DsMu",
    "DsTau",
    "Hb_bs_fd",
    "Hb_bs_dc",

    "Hb_bpm_fd",
    "Hb_bpm_dc",

    "Hb_b0_fd",
    "Hb_b0_dc",

    "Hb_lambdab_fd",
    "Hb_lambdab_dc",

    "Hb_others"]:

    #for key in ["comb","DsStarMu","DsMu","Hb","DsStarTau","DsTau"]:
    h = histos[key].Clone()
    h.Divide(hErr)
    norm_stack.Add(h)

  norm_stack.Draw('hist same')

  #subplot
  line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
  line.SetLineColor(ROOT.kBlack)
  line.SetLineWidth(1)
  ratio_stats.GetYaxis().CenterTitle()

  if logy:
    ROOT.gPad.SetLogy()

  if logx:
    ROOT.gPad.SetLogx()


  ratio_stats.Draw('E2')
  norm_stack.Draw('hist same')
  ratio_stats.Draw('E2 same')
  line.Draw('same')

  
  if dataBlind and not blind:
    pass 
  else:
    ratio.Draw('EP same')

  ROOT.gPad.RedrawAxis()
  
  #saving
  c1.Modified()
  c1.Update()


  name = ""
  folder = "/"

  if region:
    name += f"{region}" 

  if blind: 
    toSave = toSave_plots + "blind/"
  else: 
    toSave = toSave_plots

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
  if not os.path.exists(toSave + "/logx"): 
    os.makedirs(toSave + "/logx")  
  if not os.path.exists(toSave + "/logy"): 
    os.makedirs(toSave + "/logy")  

  # save nn info
  with open( toSave + f"/info.txt", "a") as f:
    f.write(f" ------ Stacked plot for variable: {var} --------- \n")
    #if pastNN: 
    #  f.write( f" These plots use the following NN model: {nnModel} \n")
    #  f.write( f" And the following score cut: {score_cut} \n")

    f.write(f"writing into info file with blind: {blind} \n")
    if not blind and not logx and not logy:
      print("Writing to:", toSave + "/info.txt")
      for key in histos.keys():
        f.write( f" Overall number of {key} entries: {histos[key].Integral()} \n")
  
        #n_bins = histos[key].GetNbinsX() 
        #for i in range(1, n_bins + 1):  # ROOT bins start from 1
        #  f.write( f" number of {key} entries in bin {i}: {histos[key].GetBinContent(i)} \n")


  if logx:
    c1.SaveAs( toSave + f"/logx/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/logx/{var}_{name}.png")

  elif logy:
    c1.SaveAs( toSave + f"/logy/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/logy/{var}_{name}.png")


  else:
    c1.SaveAs( toSave + f"/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")


  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()


  # last but not least calculate the weights between sf and data
  # this is only relevant for the first ds+mu pt stacked histo
  weights = []
  n_bins = hComb.GetNbinsX() 
  for i in range(1, n_bins + 1):  # ROOT bins start from 1

    bin_low_edge  = hComb.GetBinLowEdge(i)
    bin_high_edge = hComb.GetBinLowEdge(i + 1)

    bin_comb = hComb.GetBinContent(i)
    bin_data = histos["Data"].GetBinContent(i)
    if bin_comb > 0:
      weight   = bin_data / bin_comb
    else:
      weight = 1
    sub_str  = f" (({var} >= {bin_low_edge}) && ({var} < {bin_high_edge} )) * {weight}" 
    weights.append(sub_str)

  weight_str = " + ".join(weights)
  #print(weight_str)

  return returnHisto, weight_str



def combSys(histos2, var, scale_hb, scale_pimu = None, scale_rest = None, logx = False, logy = False, region = None, blind = False):

  histos = {}
  for key, hist in histos2.items():

    #print(key) 
    # we skip the blind keys and rename the DsTau_blind to DsTau for the corresponding cases
    if "blind" in key: continue

    if   (key == "DsTau"     and blind == True) : histos[key] = histos2["DsTau_blind"].Clone()
    elif (key == "DsStarTau" and blind == True) : histos[key] = histos2["DsStarTau_blind"].Clone()
    else                                        : histos[key] = hist.Clone()
    #print(histos) 

  color, labels = getColorAndLabelSignalDistinction("stacked")

  #supress data plotting in SR. region is a string of type: ch0
  dataBlind = False
  if ((int(region[2:]) > 6) and (int(region[2:]) < 11)): 
    dataBlind = True 
    color["DsTau"] = color["DsMu"]
    color["DsStarTau"] = color["DsStarMu"]

  # if we are anyway using the blinded histos (with a blind scale, no need do hide data)
  if (blind == True): dataBlind = False

  for i,key in enumerate(histos.keys()):

    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillColor(color[key])
    histos[key].SetLineColor(color[key])
    if key == "Data":
       histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["DsMu"].GetNbinsX()

  # error histo
  hErr = histos["DsMu"].Clone()
  hErr.Reset()
  hErr.SetName("err")
  hErr.GetXaxis().SetTitleOffset(1.3)
  hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())
  hErr.SetLineWidth(0)
  hErr.SetFillColor(ROOT.kBlack)
  hErr.SetFillStyle(3144)
  hErr.SetMarkerStyle(0)
  ###################################
  ## scale Hb to signals           ##
  ###################################
  for key in histos.keys():
    if ("Hb" in key) :
      #print(f"===> scale {key} to signals")
      histos[key].Scale(scale_hb)
  ###################################
  ## Scale Comb to data            ##
  ###################################
  # normalize
  hComb = histos["Data_sf_pimu"].Clone()
  hComb.Scale(scale_pimu)
  hComb.SetLineColor(ROOT.kGray + 1)
  hComb.SetFillColor(ROOT.kGray + 1)
  nComb = hComb.Integral()
  histos["comb"] = hComb # append to create datacard later
  ###################################
  ## Scale MC to data              ##
  ###################################
  for key in histos.keys():
    #normalize and scale the signals and hb to data
    if ('comb' not in key) and ('Data' not in key):
      #print(f"===> sclaing {key} to data")
      histos[key].Scale(scale_rest) # NEW
  
  ################################
  ## Build necessary histograms ##
  ################################
  # right sign is just data
  hRest = histos["DsMu"].Clone()
  hRest.Add(histos["DsTau"])
  hRest.Add(histos["DsStarMu"])
  hRest.Add(histos["DsStarTau"])
  hRest.Add(histos["Hb"])
  # and from data we now subtract the expected signal + bkg
  hComb_expected = histos["Data"].Clone()
  hComb_expected.Add(hRest, -1.0)
  ################################
  # Normalize to get shapes      #
  ################################
  hRest         .Scale(1.0/hRest.Integral())
  hComb         .Scale(1.0/hComb.Integral())
  hComb_expected.Scale(1.0/hComb_expected.Integral())
  hRest         .SetFillStyle(0)
  hComb         .SetFillStyle(0)
  hComb_expected.SetFillStyle(0)
  hRest         .SetLineColor(ROOT.kAzure)
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
  leg = ROOT.TLegend(.3,.70,.88,0.88)
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(42)
  leg.SetTextSize(0.035)
  leg.SetNColumns(1)
  leg.AddEntry(hRest                            ,'Signals + Hb' ,'L' )
  leg.AddEntry(hComb_expected                   ,'Comb. + Fakes Expectation (Exp.)' ,'L' )
  leg.AddEntry(hComb                            ,'Comb. + Fakes Estimation (Est.)'  ,'L' )
  #plot mainpad
  if logy:
    ROOT.gPad.SetLogy()
    hRest.GetYaxis().SetTitle('Events')
    hRest.GetYaxis().SetRangeUser(1e-3, hRest.GetMaximum()*1000)
    hRest.SetMinimum(100.001) # avoid idsplaying tons of comb bkg
  elif logx:
    ROOT.gPad.SetLogx()
    yAxisTitle = "events"
    hRest.GetYaxis().SetTitle(yAxisTitle)
    hRest.GetYaxis().SetRangeUser(1e-3, hRest.GetBinContent(hRest.GetMaximumBin())*1.5)
  else:
    maxi = max([hRest.GetBinContent(hRest.GetMaximumBin()), hComb.GetBinContent(hComb.GetMaximumBin())])
    yAxisTitle = "events"
    hRest.GetYaxis().SetTitle(yAxisTitle)
    hRest.GetYaxis().SetRangeUser(1e-3, maxi*1.8)
  hRest         .Draw("HIST")
  hComb         .Draw("HIST SAME")
  hComb_expected.Draw("HIST SAME")
  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Private Work', lumi_13TeV = f'')
  #plot ratiopad
  ratio_pad.cd()
  #ROOT.gPad.SetLogy()
  ratio = hComb.Clone()
  ratio.SetName('ratio')
  ratio.Divide(hComb_expected)
  #ratio_stats = hErr.Clone()
  #ratio_stats.SetName('ratiostats')
  #ratio_stats.Divide(hErr)
  ratio.SetMaximum(1.299) # avoid displaying 2, that overlaps with 0 in the main_pad
  ratio.SetMinimum(0.7) # and this is for symmetry
  ratio.GetYaxis().SetTitle('Est./ Exp.')
  ratio.GetXaxis().SetTitle(ratio.GetXaxis().GetTitle())
  ratio.GetYaxis().SetTitleOffset(0.5)
  ratio.GetYaxis().SetNdivisions(405)
  ratio.GetYaxis().SetTitleOffset(0.5)
  ratio.GetXaxis().SetLabelSize(3.* ratio.GetXaxis().GetLabelSize())
  ratio.GetYaxis().SetLabelSize(3.* ratio.GetYaxis().GetLabelSize())
  ratio.GetXaxis().SetTitleSize(3.* ratio.GetXaxis().GetTitleSize())
  ratio.GetYaxis().SetTitleSize(3.* ratio.GetYaxis().GetTitleSize())
  
  #subplot
  #line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
  #line.SetLineColor(ROOT.kBlack)
  #line.SetLineWidth(1)
  ratio.GetYaxis().CenterTitle()
  if logy:
    ROOT.gPad.SetLogy()
  if logx:
    ROOT.gPad.SetLogx()
  ratio      .Draw('E')
  #ratio_stats.Draw('E SAME')
  #line.Draw('same')
  ROOT.gPad.RedrawAxis()
  #saving
  c1.Modified()
  c1.Update()
  name = ""
  folder = "/"
  if region:
    name += f"{region}"
  if blind:
    toSave = toSave_plots + "blind/"
  else:
    toSave = toSave_plots + "combSys/"
  if not os.path.exists(toSave):
    os.makedirs(toSave)
  if not os.path.exists(toSave + "/logx"):
    os.makedirs(toSave + "/logx")
  if not os.path.exists(toSave + "/logy"):
    os.makedirs(toSave + "/logy")
  # save nn info
  with open( toSave + f"/info.txt", "a") as f:
    f.write(f" ------ Stacked plot for variable: {var} --------- \n")
    #if pastNN: 
    #  f.write( f" These plots use the following NN model: {nnModel} \n")
    #  f.write( f" And the following score cut: {score_cut} \n")
    f.write(f"writing into info file with blind: {blind} \n")
    if not blind and not logx and not logy:
      print("Writing to:", toSave + "/info.txt")
      for key in histos.keys():
        f.write( f" Overall number of {key} entries: {histos[key].Integral()} \n")
        #n_bins = histos[key].GetNbinsX() 
        #for i in range(1, n_bins + 1):  # ROOT bins start from 1
        #  f.write( f" number of {key} entries in bin {i}: {histos[key].GetBinContent(i)} \n")
  if logx:
    c1.SaveAs( toSave + f"/logx/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/logx/{var}_{name}.png")
  elif logy:
    c1.SaveAs( toSave + f"/logy/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/logy/{var}_{name}.png")
  else:
    c1.SaveAs( toSave + f"/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")
  
  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()
  # last but not least calculate the weights between sf and data
  # this is only relevant for the first ds+mu pt stacked histo
  weights = []
  n_bins = hComb.GetNbinsX()
  for i in range(1, n_bins + 1):  # ROOT bins start from 1
    bin_low_edge  = hComb.GetBinLowEdge(i)
    bin_high_edge = hComb.GetBinLowEdge(i + 1)
    bin_comb = hComb.GetBinContent(i)
    bin_data = histos["Data"].GetBinContent(i)
    if bin_comb > 0:
      weight   = bin_data / bin_comb
    else:
      weight = 1
    sub_str  = f" (({var} >= {bin_low_edge}) && ({var} < {bin_high_edge} )) * {weight}"
    weights.append(sub_str)
  weight_str = " + ".join(weights)

  #print(weight_str)
  return returnHisto



def shapesPlot(histos2, var, hb_scale, abc = None, scale_bkg = None, scale_kk = None, scale_pimu = None, scale_rest = None, scale_n = None, log = False, fakemass = None, A = None, B = None, C = None, S = None, bs = None, b0 = None, bplus = None, region = None, blind = False):


  histos = {}
  for key, hist in histos2.items():

    #print(key) 
    # we skip the blind keys and rename the DsTau_blind to DsTau for the corresponding cases
    if "blind" in key: continue

    if   (key == "DsTau"     and blind == True) : histos[key] = histos2["DsTau_blind"].Clone()
    elif (key == "DsStarTau" and blind == True) : histos[key] = histos2["DsStarTau_blind"].Clone()
    else                                        : histos[key] = hist.Clone()
    #print(histos) 

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    #print("my key in histos is: ", key)
    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(color[key])
    histos[key].SetLineWidth(3)

    #if key == "Data":
    #   histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["DsMu"].GetNbinsX()

  # stacked histo
  hs = ROOT.THStack(var,"")

  # error histo
  hErr = histos["DsMu"].Clone()  
  hErr.Reset()
  hErr.SetName("err")
  hErr.GetXaxis().SetTitleOffset(1.3)
  hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())  
  hErr.SetLineWidth(0)
  hErr.SetFillColor(ROOT.kBlack)
  hErr.SetFillStyle(3144)
  hErr.SetMarkerStyle(0)

  ###################################
  ## scale Hb to signals           ##
  ###################################

  # for key in histos.keys():
  #   if ("Hb" in key) : 
  #     print(f"===> scale {key} to signals")
  #     histos[key].Scale(scale_hb)

  ###################################
  ## Scale Comb to data            ##
  ###################################

  # normalize
 
  hComb = histos["Data_sf_pimu"].Clone()
  #hComb.Scale(scale_pimu)

  hComb.SetLineColor(ROOT.kGray + 1)
  hComb.SetFillColor(ROOT.kGray + 1)

  nComb = hComb.Integral()
  histos["comb"] = hComb # append to create datacard later

  ###################################
  ## Scale MC to data              ##
  ###################################

  #for key in histos.keys():

  #  #normalize and scale the signals and hb to data
  #  if ('comb' not in key) and ('Data' not in key):
  #    print(f"===> sclaing {key} to data")
  #    histos[key].Scale(scale_rest) # NEW

  ###################################
  ## Normalize                     ## 
  ###################################

  for key in histos.keys():
    try: histos[key].Scale( 1 / histos[key].Integral())
    except: pass #in case no mc (can happen if we test on small subsample)

  ###################################
  ## Create Pads                   ## 
  ###################################

  c1 = ROOT.TCanvas(var, '', 700, 700)
  c1.Draw()
  c1.cd()
  main_pad = ROOT.TPad('main_pad', '', 0., 0.0, 1. , 1.  )
  main_pad.Draw()
  c1.cd()
  #ratio_pad = ROOT.TPad('ratio_pad', '', 0., 0., 1., 0.25)
  #ratio_pad.Draw()
  main_pad.SetTicks(True)
  #main_pad.SetBottomMargin(0.)
  #main_pad.SetLeftMargin(.16)
  #ratio_pad.SetTopMargin(0.)
  #ratio_pad.SetLeftMargin(.16)
  #ratio_pad.SetGridy()
  #ratio_pad.SetBottomMargin(0.45)
  if log:
    ROOT.gPad.SetLogx()

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
  leg.SetNColumns(2)
  leg.AddEntry(histos["DsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'L' )
  leg.AddEntry(histos["DsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'L' )
  leg.AddEntry(histos["DsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','L' )
  leg.AddEntry(histos["DsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'L' )

  leg.AddEntry(histos["Hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'L' )
  leg.AddEntry(histos["Data"]       ,'Data','LEP')
  leg.AddEntry(histos["comb"]                           ,'Comb. + Fakes'  ,'L' )
  
  #plot mainpad
  key1 = next(iter(histos)) # get first key (will be drawn first)

  if log:
    ROOT.gPad.SetLogy()

  # draw all
  for i,key in enumerate(["DsTau","DsStarTau","Hb","DsMu","DsStarMu","comb","Data"]):

    #print(f"drawing variable {key}")

    if i == 0:
      histos[key].GetYaxis().SetTitle('a.u.')
      histos[key].Draw("HIST")
      if log:
        histos[key].GetYaxis().SetRangeUser(1e-3, getYMax(histos[key]) * 1000 ) 
      else:
        histos[key].GetYaxis().SetRangeUser(1e-3, getYMax(histos[key]) * 1.5 ) 
    else:
      histos[key].Draw("HIST SAME")

  leg.Draw("SAME")

  #CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  
  #saving
  c1.Modified()
  c1.Update()

  name = ""
  folder = "/"

  if region:
    name += f"{region}" 

  toSave = toSave_plots + "/shapes/"

  if blind: 
    toSave += "blind/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  if log:
    c1.SaveAs( toSave + f"/log/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/log/{var}_{name}.png")
  else:
    c1.SaveAs( toSave + f"/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")


  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()

  return returnHisto

def shapesBkg(histos2, var , region = None):

  histos = {key: hist.Clone() for key, hist in histos2.items()}
  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(color[key])
    histos[key].SetLineWidth(3)

    if key == "Data":
       histos[key].SetMarkerStyle(8) 


  # mass histo for fakemass
  histos["Data_sf_kk"]  .SetLineColor(ROOT.kGreen)
  histos["Data_sf_pimu"].SetLineColor(ROOT.kCyan)
  #histos["data_sf_both"].SetLineColor(ROOT.kBlue)
  #histos["data_sf_one"].SetLineColor(ROOT.kOrange)

  ###################################
  ## Create Pads                   ## 
  ###################################

  c1 = ROOT.TCanvas(var, '', 700, 700)
  c1.Draw()
  c1.cd()
  main_pad = ROOT.TPad('main_pad', '', 0., 0.0, 1. , 1.  )
  main_pad.Draw()
  c1.cd()
  main_pad.SetTicks(True)

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
  leg.SetNColumns(2)

  leg.AddEntry(histos["Data"]           ,'Right sign data'  ,'L' )
  leg.AddEntry(histos["Data_sf_kk"]     ,'kk flip'          ,'L' )
  leg.AddEntry(histos["Data_sf_pimu"]   ,'pimu flip'        ,'L' )
  #leg.AddEntry(histos["data_sf_both"]   ,'both AND'         ,'L' )
  #leg.AddEntry(histos["data_sf_one"]    ,'both OR'          ,'L' )
  
  #plot mainpad

  #normalize -> drawing with NORM doesnt work!!!
  histos["Data"]        .Scale(1./histos["Data"]        .Integral())  
  histos["Data_sf_kk"]  .Scale(1./histos["Data_sf_kk"]  .Integral())   
  histos["Data_sf_pimu"].Scale(1./histos["Data_sf_pimu"].Integral())  
  #histos["data_sf_both"].Scale(1./histos["data_sf_both"].Integral())
  #histos["data_sf_one"].Scale(1./histos["data_sf_one"].Integral())

  histos["Data"].GetYaxis().SetTitle('a.u.')
  histos["Data"].GetYaxis().SetRangeUser(1e-3, histos["Data"].GetMaximum()*2.0)
 
  histos["Data"]        .Draw("EP")
  histos["Data_sf_kk"]  .Draw("HIST SAME") 
  histos["Data_sf_pimu"].Draw("HIST SAME")
  #histos["data_sf_both"].Draw("HIST SAME")
  #histos["data_sf_one"].Draw("HIST SAME")

  leg.Draw("SAME")

  #saving
  c1.Modified()
  c1.Update()

  name = ""
  folder = "/"

  if region:
    name += f"{region}" 

  toSave = toSave_plots + "/shapes_bkg/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  with open( toSave + f"/info.txt", "a") as f:
    # save nn info
    #if pastNN: 
    #  f.write( f" These plots use the following NN model: {nnModel} \n")
    #  f.write( f" And the following score cut: {score_cut} \n")
    # write selection info
    # save the bin entry numbers in txt file
    for key in histos.keys():
      f.write( f" Overall number of {key} entries: {histos[key].Integral()} \n")

      n_bins = histos[key].GetNbinsX() 
      for i in range(1, n_bins + 1):  # ROOT bins start from 1
        f.write( f" number of {key} entries in bin {i}: {histos[key].GetBinContent(i)} \n")
      

  c1.SaveAs( toSave + f"/{var}_{name}.pdf")
  c1.SaveAs( toSave + f"/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")


  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()

  return returnHisto







def normPlot(histos, var, scale_kk, fakemass = None , A = None ,B = None ,C = None ,S = None ,log = False):
 
  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):

    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(color[key])
    histos[key].SetLineWidth(3)

  # take as model
  bins  = histos["DsMu"].GetNbinsX()
  start = histos["DsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  stop  = histos["DsMu"].GetXaxis().GetBinUpEdge(bins) # last bin

  # mass histo for fakemass
  hComb = ROOT.TH1D("mass","mass",bins,start,stop)
  hComb.GetXaxis().SetTitleOffset(1.3)
  hComb.GetYaxis().SetTitleSize(1*hComb.GetYaxis().GetTitleSize())
  hComb.GetXaxis().SetTitleSize(1*hComb.GetXaxis().GetTitleSize())  
  hComb.SetName("mass")
  hComb.SetFillStyle(0)
  hComb.SetLineColor(ROOT.kGray + 2)
  hComb.SetLineWidth(3)

  ###################################
  ## Combinatorial treatment       ##
  ###################################

  #simple: take sign slipped data as comb
  hComb = histos["Data_sf_pimu"].Clone()
  hComb.Scale(scale_pimu)


  try:  hComb.Scale(1 / hComb.Integral())
  except: pass #in case no comb (never happens)
  ###################################
  ## Normalize                     ##
  ###################################

  for key in histos.keys():
    try: histos[key].Scale( 1 / histos[key].Integral())
    except: pass #in case no mc (can happen if we test on small subsample)
 
  ###################################
  ## Create Pads                   ## 
  ###################################

  c1 = ROOT.TCanvas(var, '', 700, 700)
  c1.Draw()
  c1.cd()
  main_pad = ROOT.TPad('main_pad', '', 0., 0., 1. , 1.  )
  main_pad.Draw()
  c1.cd()
  main_pad.SetTicks(True)
  #main_pad.SetBottomMargin(0.)
  #main_pad.SetLeftMargin(.16)
  #ratio_pad.SetTopMargin(0.)
  #ratio_pad.SetLeftMargin(.16)
  #ratio_pad.SetGridy()
  #ratio_pad.SetBottomMargin(0.45)
  if log:
    ROOT.gPad.SetLogy()

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
  leg.SetTextSize(0.03)
  leg.SetNColumns(3)
  leg.AddEntry(histos["DsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'L' )
  leg.AddEntry(histos["DsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'L' )
  leg.AddEntry(histos["DsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','L' )
  leg.AddEntry(histos["DsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'L' )
  leg.AddEntry(histos["Hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  leg.AddEntry(hComb                           ,'Comb. + Fakes'  ,'L' )
  
  #plot mainpad
  key1 = next(iter(histos)) # get first key (will be drawn first)

  if log:
    ROOT.gPad.SetLogy()
    histos[key1].GetYaxis().SetRangeUser(1e-3, getYMax(list(histos.values()) + [hComb]) * 1000)#histos[key1].GetMaximum()*1000)

  else:
    histos[key1].GetYaxis().SetRangeUser(1e-3, getYMax(list(histos.values()) + [hComb]) * 1.5 ) #histos[key1].GetMaximum()*1.5)
    histos[key1].GetYaxis().SetTitle('a.u.')

  # draw all
  for i,key in enumerate(histos.keys()):

    print(f"drawing variable {key}")

    if ("comb" in key): 
      #comes later, dont draw sidebands
      continue

    if i == 0:
      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")

  hComb.Draw("HIST SAME") 
  leg.Draw("SAME")
  #ROOT.gPad.RedrawAxis()
  
  if trigger == "mu7": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  if trigger == "mu9": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = '20.7 fb^{-1}')

  #saving
  c1.Modified()
  c1.Update()


  name = ""
  folder = "/"

  toSave = f"/work/pahwagne/RDsTools/plots/normplots/{dt}/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  # save nn info

  if log:
    c1.SaveAs( toSave + f"/log/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/log/{var}_{name}.png")
  else:
    c1.SaveAs( toSave + f"/{var}_{name}.pdf")
    c1.SaveAs( toSave + f"/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")

def methodDistinction(histos,name, selection,norm = True):

  "Histos is a DICTIONARY!"

  #only makes directory if not already existing
  os.system(f"mkdir -p ./method_dist/")

  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  legend = getLegend("method_dist", nColumns = 3)
  canvas.cd()

  yMax = getYMax(list(histos.values()),norm)

  for i, hist in enumerate(histos.values()):

    var = list(histos.keys())[i]

    if norm: 
      hist.GetValue().Scale(1/hist.Integral())
      hist.GetValue().GetYaxis().SetTitle("a.u.")

    hist.SetMaximum(1.3*yMax)    
    _,label = getColorAndLabel(var)

    # draw 
    if  (i == 0) : hist.GetValue().Draw("HIST")
    else         : hist.GetValue().Draw("HIST SAME")

    legend.AddEntry(hist.GetValue(),label,"l")
    legend.Draw("SAME")

  toSave = "/work/pahwagne/RDsTools/plots/method_dist/{dt}/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)


  canvas.SaveAs( toSave + f"/{name}.pdf")  
  canvas.SaveAs( toSave + f"/{name}.png")  
  print(f"DONE")


def writeShapes(hist_dict, outputFile, binned = False, channel = "placeholder :)", var = None):

      hammer_sys = False
      bs_tau_sys = False
      comb_sys   = False

      for key in hist_dict.keys():
        if "e1Up"      in key: hammer_sys = True
        if "bsTauUp"   in key: bs_tau_sys = True
        if "combUp"    in key: comb_sys   = True


      #print(hist_dict["dsMu"])
      #print(outputFile)
      outputFile.WriteObject(hist_dict["DsMu"],       "dsMu"      + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["DsTau"],      "dsTau"     + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["DsStarMu"],   "dsStarMu"  + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["DsStarTau"],  "dsStarTau" + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb"],         "hb"          + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb_bs_fd"],       "hb_bs_fd"        + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["Hb_bs_dc"],       "hb_bs_dc"        + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb_b0_fd"],       "hb_b0_fd"        + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["Hb_b0_dc"],       "hb_b0_dc"        + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb_bpm_fd"],      "hb_bpm_fd"       + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["Hb_bpm_dc"],      "hb_bpm_dc"       + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb_lambdab_fd"],  "hb_lambdab_fd"   + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["Hb_lambdab_dc"],  "hb_lambdab_dc"   + binned * f"_ch{channel}" )

      outputFile.WriteObject(hist_dict["Hb_others"],      "hb_others"   + binned * f"_ch{channel}" )


      outputFile.WriteObject(hist_dict["Data"],       "data_obs"  + binned * f"_ch{channel}" ) #combine convention
      outputFile.WriteObject(hist_dict["comb"],       "comb"      + binned * f"_ch{channel}" )


      if hammer_sys:
        # if there are systematics, include them in the scaling
        for s in sys_scalar:
          for direc in ["Up", "Down"]:
            outputFile.WriteObject(hist_dict["DsTau_"        + s + direc],"dsTau_"        + s + scalar_model + direc + binned * f"_ch{channel}"  )
            outputFile.WriteObject(hist_dict["DsMu_"         + s + direc],"dsMu_"         + s + scalar_model + direc + binned * f"_ch{channel}"  )

        for s in sys_vector:
          for direc in ["Up", "Down"]:
            outputFile.WriteObject(hist_dict["DsStarTau_"    + s + direc],"dsStarTau_"    + s + vector_model + direc + binned * f"_ch{channel}"  )
            outputFile.WriteObject(hist_dict["DsStarMu_"     + s + direc],"dsStarMu_"     + s + vector_model + direc + binned * f"_ch{channel}"  )

      if bs_tau_sys:
        outputFile.WriteObject(hist_dict["DsTau_bsTauUp" ]    ,"dsTau_bsTauUp"     + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsMu_bsTauUp"  ]    ,"dsMu_bsTauUp"      + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsStarTau_bsTauUp" ],"dsStarTau_bsTauUp" + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsStarMu_bsTauUp"  ],"dsStarMu_bsTauUp"  + binned * f"_ch{channel}"  )

        outputFile.WriteObject(hist_dict["DsTau_bsTauDown" ]    ,"dsTau_bsTauDown"     + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsMu_bsTauDown"  ]    ,"dsMu_bsTauDown"      + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsStarTau_bsTauDown" ],"dsStarTau_bsTauDown" + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["DsStarMu_bsTauDown"  ],"dsStarMu_bsTauDown"  + binned * f"_ch{channel}"  )

      if comb_sys:
        outputFile.WriteObject(hist_dict["Data_sf_pimu_combUp"   ]    ,"comb_combSysUp"       + binned * f"_ch{channel}"  )
        outputFile.WriteObject(hist_dict["Data_sf_pimu_combDown" ]    ,"comb_combSysDown"     + binned * f"_ch{channel}"  )


      #myFile_1D.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      #myFile_1D.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      #myFile_1D.WriteObject(histosScaled["hb"],        "hb"        )
      #myFile_1D.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      #myFile_1D.WriteObject(histosScaled["comb"],      "comb"      )


def writeBinnedDatacard(histos, var, region, digits = 5, blind = False):

  name = ""


  if blind:
    toSave = toSave_cards + "/blind"
  else:
    toSave = toSave_cards

  if not os.path.exists(toSave): 
    os.makedirs(toSave)

  hammer_sys = False
  bs_tau_sys = False
  comb_sys   = False

  for key in histos.keys():
    if "e1Up" in key: hammer_sys = True
    if "bsTauUp" in key: bs_tau_sys = True
    if "combUp" in key: comb_sys = True


  if hammer_sys and not bs_tau_sys:
    #if args.findcut:
    #  temp = open("/work/pahwagne/RDsTools/fit/datacardTemplateSystematics_binned_cuts.txt", "rt")
    #temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned.txt", "rt")
    #temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_starVSnonstar_hbYields.txt", "rt")
    temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_starVSnonstar_hbYields_comb.txt", "rt")
    #temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_hb_varies.txt", "rt")
    #temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_reducedFF.txt", "rt")
    #temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_old.txt", "rt")

  elif hammer_sys and bs_tau_sys:
    temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_reducedFF_and_bsTau.txt", "rt")
  elif not hammer_sys and bs_tau_sys:
    temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplateSystematics_binned_bsTau.txt", "rt")

  else:
    #if args.findcut:
    #  temp = open("/work/pahwagne/RDsTools/fit/datacardTemplate_binned_cuts.txt", "rt")
    temp = open("/work/pahwagne/RDsTools/fit/new_temps/datacardTemplate_binned.txt", "rt")

  card_dir = toSave + f"/datacard_binned_{var}_ch{region}.txt"
  card = open(card_dir, "wt")

  # width in datacard template is 17 spaces
  spaces = 22

  dataStr       = str(round(histos["Data"].Integral(), digits))
  dataStr      += " "*(spaces - len(dataStr))

  dsMuStr       = str(round(histos["DsMu"].Integral(), digits))
  dsMuStr      += " "*(spaces - len(dsMuStr))

  dsTauStr      = str(round(histos["DsTau"].Integral(), digits))
  dsTauStr     += " "*(spaces - len(dsTauStr))

  dsStarMuStr   = str(round(histos["DsStarMu"].Integral(), digits))
  dsStarMuStr  += " "*(spaces - len(dsStarMuStr))

  dsStarTauStr  = str(round(histos["DsStarTau"].Integral(), digits))
  dsStarTauStr += " "*(spaces - len(dsStarTauStr))

  #hbStr         = str(round(histos["Hb"].Integral(), digits))
  #hbStr        += " "*(spaces - len(hbStr))

  hb_bs_fd_Str  = str(round(histos["Hb_bs_fd"].Integral(), digits))
  hb_bs_fd_Str += " "*(spaces - len(hb_bs_fd_Str))

  hb_bs_dc_Str  = str(round(histos["Hb_bs_dc"].Integral(), digits))
  hb_bs_dc_Str += " "*(spaces - len(hb_bs_dc_Str))

  hb_bpm_fd_Str  = str(round(histos["Hb_bpm_fd"].Integral(), digits))
  hb_bpm_fd_Str += " "*(spaces - len(hb_bpm_fd_Str))

  hb_bpm_dc_Str  = str(round(histos["Hb_bpm_dc"].Integral(), digits))
  hb_bpm_dc_Str += " "*(spaces - len(hb_bpm_dc_Str))

  hb_b0_fd_Str  = str(round(histos["Hb_b0_fd"].Integral(), digits))
  hb_b0_fd_Str += " "*(spaces - len(hb_b0_fd_Str))

  hb_b0_dc_Str  = str(round(histos["Hb_b0_dc"].Integral(), digits))
  hb_b0_dc_Str += " "*(spaces - len(hb_b0_dc_Str))

  hb_lambdab_fd_Str  = str(round(histos["Hb_lambdab_fd"].Integral(), digits))
  hb_lambdab_fd_Str += " "*(spaces - len(hb_lambdab_fd_Str))

  hb_lambdab_dc_Str  = str(round(histos["Hb_lambdab_dc"].Integral(), digits))
  hb_lambdab_dc_Str += " "*(spaces - len(hb_lambdab_dc_Str))

  hb_others_Str  = str(round(histos["Hb_others"].Integral(), digits))
  hb_others_Str += " "*(spaces - len(hb_others_Str))



  combStr       = str(round(histos["comb"].Integral(), digits))  
  #combStr      += " "*(spaces - len(combStr))

  #rates = dsMuStr + dsTauStr + dsStarMuStr + dsStarTauStr + hbStr + combStr
  rates = dsMuStr  + \
          dsTauStr + \
          dsStarMuStr + \
          dsStarTauStr + \
          hb_bs_fd_Str + \
          hb_bs_dc_Str + \
          hb_bpm_fd_Str + \
          hb_bpm_dc_Str + \
          hb_b0_fd_Str + \
          hb_b0_dc_Str + \
          hb_lambdab_fd_Str + \
          hb_lambdab_dc_Str + \
          hb_others_Str + \
          combStr

  for line in temp:
    if "HOOK_REGION"            in line: line = line.replace("HOOK_REGION",    str(region))
    if "HOOK_RATES"             in line: line = line.replace("HOOK_RATES",     rates      )
    if "HOOK_DATA_RATE"         in line: line = line.replace("HOOK_DATA_RATE", dataStr    )
    if "HOOK_VAR"               in line: line = line.replace("HOOK_VAR",       var        )
    #if "HOOK_NAME"              in line: line = line.replace("HOOK_NAME",      name       )
    if "HOOK_DATETIME"          in line: line = line.replace("HOOK_DATETIME",  dt )
    if "HOOK_CUT"               in line: line = line.replace("HOOK_CUT",       args.cut )
    if "HOOK_BLIND"             in line: line = line.replace("HOOK_BLIND",     blind * "blind" )
    #if "HOOK_TRIGGER"           in line: print("IM HERE!"); line = line.replace("HOOK_TRIGGER",   trigger  )
    card.write(line)

  temp.close()
  card.close()

  #return the path to the card
  return card_dir



def loadHistos(path,channel,i):

  f = ROOT.TFile.Open(f"{path}/histos_{channel}_{i}.root")
  
  hists = {}
  #loop over all keys in the root file
  for key in f.GetListOfKeys():
    name = key.GetName()
    obj = f.Get(name)

    #only consider TH1, TH1F   
    if not obj or not obj.InheritsFrom("TH1"): continue

    obj.SetDirectory(0)
   
    hists[name] = obj

  f.Close()

  return hists



###########################################
# Import histograms, one for every region #
###########################################

# import the histograms from the root file!

##########################
# Extract number of bins #
##########################

# root filenames: histos_channel_X.root, where X is the bin. 
# take channel DsMu as example to extract the # of bins

# logic \d means any digit (0-9)
# \d+ means any higher number (including 10,11, ...)
# () the pattern is always placed into parantheses
# \. means "dot" in regex, otherwise only . means "any character" 
pattern = re.compile(r"histos_DsMu_(\d+)\.root")
bins = []

#loop over all files
for f in os.listdir(toSave_plots):
  #if file matches the pattern extract the pattern (group(1))
  match = pattern.match(f)
  if match:
    bins.append(int(match.group(1)))
  
#get the maximum 
num_bins = max(bins)

############################
# Start plotting loop      #
############################
#load from json file
with open(f"{toSave_plots}/normalization_parameters.json", "r") as f:
  data = json.load(f)
  print(f" ====> JSON file loaded: {data}")
  hb_ratio          = data["hb_ratio_massfit"]
  hb_scale_massfit  = data["hb_scale_massfit"]
  scale_pimu        = data["scale_pimu"]
  scale_pimu_blind  = data["scale_pimu_blind"]
  scale_rest        = data["scale_rest"]
  scale_rest_blind  = data["scale_rest_blind"]
  #trigger           = data["trigger"]
 
#for every variable we plot, we get the total number of Hb and Mu events to get the Hb scale

#events_Hb            = {}
#events_Mu_in_Hb      = {}
#events_DsMu_woHammer = {}
#
#for var in vars_base:
#
#  #create dict entry
#  events_Hb           [var] = 0
#  events_DsMu_woHammer[var] = 0
#  events_Mu_in_Hb     [var] = 0
#
#  for i in range(num_bins+1):
#
#    #load Hb and DsMu histo for this channel
#    chan_histos_Hb       = loadHistos(toSave_plots, "Hb",   i)
#    chan_histos_Mu_in_Hb = loadHistos(toSave_plots, "Mu_in_Hb", i)
#    chan_histos_DsMu_woHammer     = loadHistos(toSave_plots, "DsMu_woHammer", i)
#
#    #load histo for this channel and var
#    events_Hb           [var] += chan_histos_Hb           [var].Integral()
#    events_Mu_in_Hb     [var] += chan_histos_Mu_in_Hb     [var].Integral()
#    events_DsMu_woHammer[var] += chan_histos_DsMu_woHammer[var].Integral()



for i in range(num_bins+1):

  writeThisBin = True

  ###############################
  # Extract variable base names #
  ###############################
  
  #pattern to look for
  #pattern = re.compile(r"^(?P<var>.+?)(?:_e\d+(?:Up|Down))?$")
  pattern = re.compile(r"^(?P<var>.+?)(?:_[^_]+(?:Up|Down))?$") 
 
  #create a set to automatically reduce
  vars_base = set()
  
  #since we can have diff variables for different channels,
  #we need to compute for every channel explicitly!
  
  # take channel DsMu as example
  fstr = f"{toSave_plots}/histos_DsMu_{i}.root"
  f = ROOT.TFile.Open(fstr)
  
  #contains files of type:
  # phiPi_m;1	
  # phiPi_m_e1Up;1	
  # phiPi_m_e1Down;1
  # ...
  # ---> extract only phiPi_m
  
  #loop over keys and match to pattern
  for key in f.GetListOfKeys():
      name = key.GetName()
  
      m = pattern.match(name)
      if not m: continue
      vars_base.add(m.group("var"))
  
  f.Close()


  
    


  channels = [
    "DsTau",
    "DsTau_blind",
    "DsStarTau",
    "DsStarTau_blind",
    "DsMu",
    "DsMu_woHammer",
    "DsStarMu",
    "Hb",
    "Hb_fd",
    "Hb_dc",

    "Hb_bs",
    "Hb_bs_fd",
    "Hb_bs_dc",

    "Hb_bpm",
    "Hb_bpm_fd",
    "Hb_bpm_dc",

    "Hb_b0",
    "Hb_b0_fd",
    "Hb_b0_dc",

    "Hb_lambdab",
    "Hb_lambdab_fd",
    "Hb_lambdab_dc",

    "Hb_others",
    "Data",
    "Data_sf_pimu",
  ]

  histos = {}

  #loop over channels
  for chan in channels:

    #contains all variables (including Up/Down) for this channel
    chan_histos = loadHistos(toSave_plots, chan, i)

    #add to big histogram
    histos[chan] = chan_histos 


  ###########################
  # Reshuffle naming logic  #
  ###########################

  from collections import defaultdict
  histos_reshuffled = defaultdict(dict)

  #remove 
  channels.remove("DsMu_woHammer")
  
  #loop over all channels
  for chan in channels:
    # and all variable bases: phiPi_m, q2_coll, ...
    for var in vars_base:
      #and the keys of histos[chan]
      for key in histos[chan].keys():
        if var in key: #meaning key must be like phiPi_m, phiPi_m_e1Up, ... etc.
          #if "Up" in key or "Down" in key:
          #else: 
          pattern = re.compile(rf'^({re.escape(var)})(?:_(.+(?:Up|Down).*))?$')
          match = pattern.match(key)
          if not match:
            #happens f.e. when var is lxy_ds and key is lxy_ds_sig
            continue

          if match.group(2):
            #Up/Down found
            histos_reshuffled[var][f"{chan}_{match.group(2)}"] = histos[chan][key]
          else: 
            histos_reshuffled[var][f"{chan}"] = histos[chan][key]


  scale_hb = hb_scale_massfit

  cards       = {var: [] for var in vars_base}
  cards_blind = {var: [] for var in vars_base}

  for var in vars_base:

    hammer_sys = False
    for key in histos_reshuffled[var].keys():
      if "e1Up"    in key: hammer_sys = True
  
  
    if var == "phiPi_m" and i == 0 and hammer_sys:
  
      answer = input("Skip ch0 and use it for yield estimation only (mass hammer systematics)? (y/n): ").strip().lower()
      if answer == "y":
        print("You chose yes!")
  
        print("##################################################################") 
        print("# For the mass, we are freezing the shape and only change yields #") 
        print("##################################################################") 
  
        modifyHammer = True
        writeThisBin = False #ch0 here has only the purpose of re-defining hammer shapes, no datacard needed

        #save yields globally
        for s in sys_scalar:
          for direc in ["Up", "Down"]:
            massYields["DsTau_"    + s + direc] = histos_reshuffled[var]["DsTau_"    + s + direc].Integral()/ histos_reshuffled[var]["DsTau"].Integral() 
            massYields["DsMu_"     + s + direc] = histos_reshuffled[var]["DsMu_"     + s + direc].Integral()/ histos_reshuffled[var]["DsMu" ].Integral() 
                                      
        for s in sys_vector:          
          for direc in ["Up", "Down"]:
            massYields["DsStarTau_"+ s + direc] = histos_reshuffled[var]["DsStarTau_"+ s + direc].Integral()/ histos_reshuffled[var]["DsStarTau"].Integral() 
            massYields["DsStarMu_" + s + direc] = histos_reshuffled[var]["DsStarMu_" + s + direc].Integral()/ histos_reshuffled[var]["DsStarMu" ].Integral() 
        print(massYields)
      
      elif answer == "n":
        print("You chose no!")
      else:
        print("Invalid input!")
   
 
  
    if var == "phiPi_m" and i == 1 and hammer_sys and modifyHammer:
 
        print("====> Redefining the ch1 (mass) hammer sys uncert.)") 
        for s in sys_scalar:
          for direc in ["Up", "Down"]:
  
            clone = histos_reshuffled[var]["DsTau"    ].Clone()
            clone.Scale( massYields["DsTau_"    + s + direc]) 
            histos_reshuffled[var]["DsTau_"    + s + direc] = clone 
  
            clone = histos_reshuffled[var]["DsMu"     ].Clone()
            clone.Scale( massYields["DsMu_"     + s + direc])
            histos_reshuffled[var]["DsMu_"     + s + direc] = clone
                                                                        
        for s in sys_vector:                                            
          for direc in ["Up", "Down"]:                                  
  
            clone = histos_reshuffled[var]["DsStarTau"    ].Clone()
            clone.Scale( massYields["DsStarTau_"+ s + direc])
            histos_reshuffled[var]["DsStarTau_"+ s + direc] = clone
  
            clone = histos_reshuffled[var]["DsStarMu"     ].Clone()
            clone.Scale( massYields["DsStarMu_" + s + direc])
            histos_reshuffled[var]["DsStarMu_" + s + direc] = clone
  
 
    #create root file for every variable which holds shapes
    shapes_name = f"/{var}_shapes_ch{i}.root"
    
    myFile       = ROOT.TFile.Open(shapes_folder       + shapes_name, "RECREATE")
    myFile_blind = ROOT.TFile.Open(shapes_folder_blind + shapes_name, "RECREATE")
    print(f"===> Producing stacked plot for variable: {var} in channel {i}") 

    histosScaled        , _ = stackedPlot(histos_reshuffled[var],       var, scale_hb,   scale_pimu = scale_pimu,       scale_rest = scale_rest,        region = f"ch{i}", blind = False )
    histosScaled_blind  , _ = stackedPlot(histos_reshuffled[var],       var, scale_hb,   scale_pimu = scale_pimu_blind, scale_rest = scale_rest_blind,  region = f"ch{i}", blind = True  )
    histosScaledLogx    , _ = stackedPlot(histos_reshuffled[var],       var, scale_hb,   scale_pimu = scale_pimu,       scale_rest = scale_rest,        region = f"ch{i}", blind = False, logx = True )
    histosScaledLogy    , _ = stackedPlot(histos_reshuffled[var],       var, scale_hb,   scale_pimu = scale_pimu,       scale_rest = scale_rest,        region = f"ch{i}", blind = False, logy = True )
  
    histosScaled_shapes     = shapesPlot(histos_reshuffled[var],        var, scale_hb,   scale_pimu = scale_pimu,       scale_rest = scale_rest,        region = f"ch{i}", blind = False )
    histosScaled_comb       = combSys   (histos_reshuffled[var],        var, scale_hb,   scale_pimu = scale_pimu,       scale_rest = scale_rest,        region = f"ch{i}", blind = False )


    ###############################################################
    # Correct one-sided hammer effects when the one-sided-ness    #
    # is less than one event (comging from numerical precision!   #
    ###############################################################


    tol = 0.001
    if hammer_sys :
 
        # get nBins wlog from central curve
        nBinsX = histosScaled["DsTau"].GetNbinsX()

        for s in sys_scalar:
          for j in range(1,nBinsX+1):

            nC    = histosScaled["DsTau"                 ].GetBinContent(j) 
            nUp   = histosScaled["DsTau_"    + s + "Up"  ].GetBinContent(j) 
            nDown = histosScaled["DsTau_"    + s + "Down"].GetBinContent(j) 

            if ( ((nC < nDown) and (nC < nUp)) or ((nC > nDown) and (nC > nUp)) ):
              print(f"===> Detected one sided error for direction {s} in bin {j} of dsTau due to numerical precision! \nFix it!") 
              if abs((nC - nUp)   / nC ) < tol: histosScaled["DsTau_"    + s + "Up"  ].SetBinContent(j,nC)
              if abs((nC - nDown) / nC ) < tol: histosScaled["DsTau_"    + s + "Down"].SetBinContent(j,nC)


            nC    = histosScaled["DsMu"                 ].GetBinContent(j) 
            nUp   = histosScaled["DsMu_"    + s + "Up"  ].GetBinContent(j) 
            nDown = histosScaled["DsMu_"    + s + "Down"].GetBinContent(j) 
                                                                                                                                   
            if ( ((nC < nDown) and (nC < nUp)) or ((nC > nDown) and (nC > nUp)) ):
              print(f"===> Detected one sided error for direction {s} in bin {j} of dsMu due to numerical precision! \nFix it!") 
              if abs((nC - nUp)   / nC ) < tol: histosScaled["DsMu_"    + s + "Up"  ].SetBinContent(j,nC)
              if abs((nC - nDown) / nC ) < tol: histosScaled["DsMu_"    + s + "Down"].SetBinContent(j,nC)


                                                                        
        for s in sys_vector:                                            
          for j in range(1,nBinsX+1):
  
            nC    = histosScaled["DsStarTau"                 ].GetBinContent(j) 
            nUp   = histosScaled["DsStarTau_"    + s + "Up"  ].GetBinContent(j) 
            nDown = histosScaled["DsStarTau_"    + s + "Down"].GetBinContent(j) 

            if ( ((nC < nDown) and (nC < nUp)) or ((nC > nDown) and (nC > nUp)) ):
              print(f"===> Detected one sided error for direction {s} in bin {j} of dsTau due to numerical precision! \nFix it!") 
              if abs((nC - nUp)   / nC ) < tol: histosScaled["DsStarTau_"    + s + "Up"  ].SetBinContent(j,nC)
              if abs((nC - nDown) / nC ) < tol: histosScaled["DsStarTau_"    + s + "Down"].SetBinContent(j,nC)



            nC    = histosScaled["DsStarMu"                 ].GetBinContent(j) 
            nUp   = histosScaled["DsStarMu_"    + s + "Up"  ].GetBinContent(j) 
            nDown = histosScaled["DsStarMu_"    + s + "Down"].GetBinContent(j) 
                                                                                                                                   
            if ( ((nC < nDown) and (nC < nUp)) or ((nC > nDown) and (nC > nUp)) ):
              print(f"===> Detected one sided error for direction {s} in bin {j} of dsMu due to numerical precision! \nFix it!") 
              if abs((nC - nUp)   / nC ) < tol: histosScaled["DsStarMu_"    + s + "Up"  ].SetBinContent(j,nC)
              if abs((nC - nDown) / nC ) < tol: histosScaled["DsStarMu_"    + s + "Down"].SetBinContent(j,nC)


 
   
    if writeThisBin:   

      print("====> Writing shapes and datacards")
      writeShapes(histosScaled,       myFile,       binned = True, channel = i, var = var)
      writeShapes(histosScaled_blind, myFile_blind, binned = True, channel = i, var = var)
    
      card       = writeBinnedDatacard(histosScaled      , var, i)
      card_blind = writeBinnedDatacard(histosScaled_blind, var, i, blind = True)
    
      cards[var].append(card)
      cards_blind[var].append(card_blind)
 
    else:

      print(f"====> Skipping region {i}, not saving any datacard ...")
 
  
  #
  ###more advanced, create binned datacards
  ##if control == "highmass":
  ##  createBinnedPlots(split,binning, controlPlotsHighMass = True) 
  ##elif control == "leftsb":
  ##  createBinnedPlots(split,binning, controlPlotsLeftSideband = True) 
  ##elif control == "rightsb":
  ##  createBinnedPlots(split,binning, controlPlotsRightSideband = True) 
  ##elif control == "complete":
  ##  createBinnedPlots(split,binning, controlPlotsComplete = True) 
  ##elif control == "custom":
  ##  createBinnedPlots(split,binning, controlPlotsCustom = True)
  ##else:
  ##  createBinnedPlots(split,binning) 
  ##
  ##
