import ROOT
import argparse
import os
import sys

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from sidebands import getSigma, getABCS
from helper import bsMass_, dsMass_, phiMass_
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

sig   = "13_05_2024_14_39_06"                      #"26_04_2024_16_28_22" #old MC 
hb    = "26_04_2024_16_28_58" #"10_04_2024_00_32_53"#"09_04_2024_18_44_36"  #"22_03_2024_17_34_53" #old inclusive
data  = "26_04_2024_18_08_15" #data -> apply SB method
bs    = "21_05_2024_20_27_35"
bplus = "21_05_2024_20_27_47"
b0    = "21_05_2024_19_58_59"
lambdab    = "21_05_2024_20_28_00"

nSignalRegion    = 3 #how many sigma until the signal region stops
nSidebands     = 5 #how many sigma until the sb starts
sbWidth = 1 # sb width

#dsMass_ = 1.96834
#bsMass_ = 5.36688
#phiMass_ = 1.019461
#########################################
## SELECTIONS                          ##
#########################################


## BASIC

selBasic        = f" (dsMu_m < {bsMass_}) & (k1_charge*k2_charge < 0) & (mu_charge*pi_charge <0) & (gen_match_success == 1))"
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

addOn1 = ' & '.join([
'(bs_pt_reco_2 > 10)',
'(cosMuW_lhcb_alt > -0.5)',
'((cosPiK1 < -0.6) || (cosPiK1 > 0.6))',
'(e_star_reco_2 < 2.5)'
])

addOn2 = ' & '.join([
'(pt_miss_lhcb_alt < 25)',
'(m2_miss_lhcb_alt > -2)',
f'(abs(kk_m - {phiMass_}) < 0.010 )',
])

addOn3 = ' & '.join([
'(m2_miss_lhcb_alt > -2)',
'(e_star_reco_weighted < 1.7)',
'(bs_pt_reco_2 > 20)',
'((cosPiK1 < -0.4) || (cosPiK1 > 0.4))',
])

bkg = ' & '.join([
f'dsMu_m > {bsMass_}',
#'mu_rel_iso_03 > 0.3'
])


#baseline = baseline #+ '&' + addOn1 + '&' + addOn2
baseline = baseline + '&' + addOn3


baselineHb =        baseline + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 5) && (gen_sig != 6) & (gen_match_success ==1 ) "
baselineDsMu =      baseline + " && (gen_sig == 0)" 
baselineDsTau =     baseline + " && (gen_sig == 1)"
baselineDsStarMu =  baseline + " && (gen_sig == 5)" 
baselineDsStarTau = baseline + " && (gen_sig == 6)"

#########################################
## CREATE RDF FROM TREE                ##
#########################################

def getRdf(dateTime, debug = False, skimmed = ""):

  if skimmed != "":
    files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_{skimmed}_{dateTime}.root" #data skimmed with selection 'skimmed'
    #files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_bkg_{dateTime}.root"  # data skimmed with kkpimu > Bs for closure
    print("taking skimmed file ...", files)
  else:
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
chainSigSB, rdfSigSB   = getRdf(sig, skimmed = "baseline"   ) # FOR SB FIT
chainSig, rdfSig   = getRdf(sig, skimmed = "baseline"   )#,  debug = True)
chainHb,  rdfHb    = getRdf(hb,  skimmed = "baseline"    )#,   debug = True)
#chainData,rdfData  = getRdf(data)#, debug = True)
chainData,rdfData  = getRdf(data,skimmed = "baseline")  #skimmed = "baseline") #pick already skimmed file!

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
    colors = {"dsMu":   ROOT.kBlue - 2,
              "dsTau":  ROOT.kGreen,
              "dsStarMu":   ROOT.kCyan,
              "dsStarTau":  ROOT.kOrange,
              "hb":         ROOT.kRed,
              "combL":      ROOT.kGray+1,
              "combR":      ROOT.kGray+1,
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
  #ntot   = chain.GetEntries( f" (gen_match_success == 1)  && ({name}_reco_1 == {name}_reco_1)")
  nReco1 = chain.GetEntries(selection + f"&& ({name}_reco_1 > {start}) && ({name}_reco_1 < {stop}) && (gen_{name} > {start}) && (gen_{name} < {stop}) && (abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && ({name}_reco_1 == {name}_reco_1)")
  #nReco1 = chain.GetEntries( f" (gen_match_success == 1) && (abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && ({name}_reco_1 == {name}_reco_1)")
    
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

def createHistos(selection,rdf, linewidth = 2, gen = True):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"

  histos = {}
   
  for var,model in models.items(): 
      print("Filling for variable ", var)
      if "gen" in var and not gen:
        #skip gen variables
        print("This is a gen variable... skip!")
        continue

      tofill = var

      print("I am filling the histo")
      if var == "phiPi_m":
        tofill = f"(run==1) * 0.9995* phiPi_m + (run!=1) * phiPi_m"
        print("correcting phiPi mass..")
        histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], "m_corr")
      else:
        histos[var] = rdf.Filter(selection).Histo1D(model[0], tofill)

      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetMaximum())
      
      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(linewidth)
  print(color)
  print(f"      Selection: {selection} DONE")

  return histos

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

  if (pureHb.Integral() == 0): 
    #avoid dividing by zero (probably muons are also empty!)
    scale_hb = 0
  else:
    scale_hb = muNew.Integral() * pureHb.Integral() / pureMu.Integral()
  print( "#Hb events for MC sample is:", scale_hb)
  return scale_hb

def stackedPlot(histos, var, hb_scale, fakemass, A,B,C,D, mlow, mhigh, log = False):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    print("adapting fill style..", key)
    histos[key].SetFillColor(color[key])
    histos[key].SetLineColor(color[key])
    if key == "data":
       histos[key].SetMarkerStyle(8) 

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
  hComb.SetLineColor(ROOT.kGray)
  # special error fillstyle 
  hErr.SetFillColor(ROOT.kGray+1)
  hErr.SetFillStyle(3344)

  #scale hb to other mc

  print("Scale of Hb: ", hb_scale )
  if hb_scale == 0: 
    histos["hb"].Scale(0.0 )
  else:
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
      hErr.Add(histos[key])
      hs.Add(histos[key])

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
  leg.AddEntry(histos["dsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'F' )
  leg.AddEntry(histos["dsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'F' )
  leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  leg.AddEntry(histos["data"]       ,'Data','LEP')
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','F' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'F' )
  leg.AddEntry(hComb                           ,'Comb. Bkg.'  ,'F' )
  leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
  #plot mainpad
  if log:
    ROOT.gPad.SetLogy()
    histos["data"].GetYaxis().SetTitle('Events')
    histos["data"].GetYaxis().SetRangeUser(1e-3, hs.GetMaximum()*1000)
    histos["data"].SetMinimum(1.001) # avoid idsplaying tons of comb bkg
  else:
    histos["data"].GetYaxis().SetRangeUser(1e-3, histos["data"].GetMaximum()*1.5)
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
  ratio_stats = histos["data"].Clone()
  ratio_stats.SetName(ratio.GetName()+'_ratiostats')
  ratio_stats.Divide(hErr)
  ratio_stats.SetMaximum(1.999) # avoid displaying 2, that overlaps with 0 in the main_pad
  ratio_stats.SetMinimum(0.0001) # and this is for symmetry
  ratio_stats.GetYaxis().SetTitle('Data / MC')
  ratio_stats.GetYaxis().SetTitleOffset(0.5)
  ratio_stats.GetYaxis().SetNdivisions(405)
  ratio_stats.GetYaxis().SetTitleOffset(0.5)
  ratio_stats.GetXaxis().SetLabelSize(3.* ratio.GetXaxis().GetLabelSize())
  ratio_stats.GetYaxis().SetLabelSize(3.* ratio.GetYaxis().GetLabelSize())
  ratio_stats.GetXaxis().SetTitleSize(3.* ratio.GetXaxis().GetTitleSize())
  ratio_stats.GetYaxis().SetTitleSize(3.* ratio.GetYaxis().GetTitleSize())
  
         
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

  print(f"I count {hComb.Integral()} comb. bkg. events for the variable {var}")
 
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


  toSave = "/work/pahwagne/RDsTools/plots/cmsplots/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/log/{var}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/log/{var}.png")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{var}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{var}.png")
  print(f"===> Produced plot: {var}.pdf")

def normPlot(histos, var, fakemass, A,B,C,S,log = False):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(color[key])
    histos[key].SetLineWidth(3)

  # take as model
  bins  = histos["dsMu"].GetNbinsX()
  start = histos["dsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  stop  = histos["dsMu"].GetXaxis().GetBinUpEdge(bins) # last bin

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

  if var != "phiPi_m":
    # use sideband method
    # left one
    hComb = histos["combL"].Clone()
    hComb.SetName("combS") #signal region
    hComb.Scale(B/(2*A))

    # right one
    hCombR = histos["combR"].Clone()
    hCombR.Scale(B/(2*C))
    hComb.Add(hCombR)

  else:
    for i in range(bins):
      hComb.SetBinContent(i+1,fakemass[i])
      hComb.SetBinError(i+1,np.sqrt(fakemass[i]))

  hComb.Scale(1 / hComb.Integral())
  ###################################
  ## Normalize                     ##
  ###################################

  for key in histos.keys():
      histos[key].Scale( 1 / histos[key].Integral())

 
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
  leg.AddEntry(histos["dsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'L' )
  leg.AddEntry(histos["dsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'L' )
  leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'L' )
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','L' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'L' )
  leg.AddEntry(hComb                           ,'Comb. Bkg.'  ,'L' )
  
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
    if "comb" in key: 
      #comes later, dont draw sidebands
      continue
    print(f"drawing variable {key}")
    if i == 0:
      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")
  hComb.Draw("HIST SAME") 
  leg.Draw("SAME")
  #ROOT.gPad.RedrawAxis()
  
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = 'X fb^{-1}')
  #saving
  c1.Modified()
  c1.Update()

  toSave = "/work/pahwagne/RDsTools/plots/normplots/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  


  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/log/{var}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/log/{var}.png")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{var}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{var}.png")
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

  toSave = "/work/pahwagne/RDsTools/plots/method_dist/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)

  canvas.SaveAs(f"/work/pahwagne/RDsTools/plots/method_dist/{name}.pdf")  
  canvas.SaveAs(f"/work/pahwagne/RDsTools/plots/method_dist/{name}.png")  
  print(f"DONE")

def createPlots(selec):

  selecHb =        selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 5) && (gen_sig != 6) & (gen_match_success ==1 ) "
  selecDsMu =      selec + " && (gen_sig == 0)" 
  selecDsTau =     selec + " && (gen_sig == 1)"
  selecDsStarMu =  selec + " && (gen_sig == 5)" 
  selecDsStarTau = selec + " && (gen_sig == 6)"

  sigma, h          = getSigma(rdfSigSB, "phiPi_m", selec + "& gen_sig == 0")
                                                                           #bins for fit #bins for fakemass
  A, B, C, S        = getABCS( data, rdfData, selec , "phiPi_m", sigma, h, binsFake = 21, nSig = nSignalRegion, nSb = nSidebands, sbWidth = sbWidth)
  #A = 24158.821588264094; B = 226543.39752377832; C = 21219.183790935822; S = 19829.186060030886; #use for debugging (faster) 
   
  #signal region
  mlow   = dsMass_ - nSignalRegion*sigma
  mhigh  = dsMass_ + nSignalRegion*sigma

  #sideband start
  mlow2  = dsMass_ - nSidebands*sigma
  mhigh2 = dsMass_ + nSidebands*sigma

  #sideband stops
  mlow3  = mlow2  - sbWidth*sigma
  mhigh3 = mhigh2 + sbWidth*sigma
  
  signalRegion = f"& ({mlow} < phiPi_m) & (phiPi_m < {mhigh})"
  leftSB       = f"& ({mlow3} < phiPi_m) & (phiPi_m < {mlow2})"
  rightSB      = f"& ({mhigh2} < phiPi_m) & (phiPi_m < {mhigh3})"
  
  # get fakemass histo
  fakemass = np.genfromtxt('mass_bincontent.csv', delimiter=',')
  
  #create histos returns a dictionary !:)
  
  #for all variables except mass plot only signal region (indicated with 'S')
  selec_S_DsMu      = createHistos(selecDsMu      + signalRegion, rdfSig , gen = False)
  selec_S_DsTau     = createHistos(selecDsTau     + signalRegion, rdfSig , gen = False)
  selec_S_DsStarMu  = createHistos(selecDsStarMu  + signalRegion, rdfSig , gen = False)
  selec_S_DsStarTau = createHistos(selecDsStarTau + signalRegion, rdfSig , gen = False)
  selec_S_Hb        = createHistos(selecHb        + signalRegion, rdfHb , gen = False)
  selec_S_Data      = createHistos(selec         + signalRegion, rdfData , gen = False)
  selecDataL        = createHistos(selec          + leftSB, rdfData , gen = False)
  selecDataR        = createHistos(selec          + rightSB, rdfData , gen = False)

  hb_scale = getHbScale(selecHb + signalRegion, selecDsMu + signalRegion)
  
  # for the ds mass plot we also want to plot the sidebands! (indicated with 'C' for complete)
  selec_C_DsMu      = createHistos(selecDsMu      , rdfSig , gen = False)
  selec_C_DsTau     = createHistos(selecDsTau     , rdfSig , gen = False)
  selec_C_DsStarMu  = createHistos(selecDsStarMu  , rdfSig , gen = False)
  selec_C_DsStarTau = createHistos(selecDsStarTau , rdfSig , gen = False)
  selec_C_Hb        = createHistos(selecHb        , rdfHb  , gen = False)
  selec_C_Data      = createHistos(selec          , rdfData , gen = False)
  
  hb_scale = getHbScale(selecHb, selecDsMu)

  for var in models.keys():
    print(f"===> Producing stacked plot for variable: {var}") 
    if "gen" in var:
      #skip gen variables
      print("This is a gen variable... skip!")
      continue
  
    if var != "phiPi_m":
  
      histos = {"dsTau":     selec_S_DsTau[var],
               "dsStarTau": selec_S_DsStarTau[var], 
               "dsMu":     selec_S_DsMu[var], 
               "dsStarMu":  selec_S_DsStarMu[var], 
               "hb"       : selec_S_Hb[var], 
               "combL"    : selecDataL[var], 
               "combR"    : selecDataR[var], 
               "data"     : selec_S_Data[var]}
      stackedPlot(histos, var, hb_scale, fakemass, A,B,C,S, mlow, mhigh)
      #stackedPlot(histos, var, hb_scale, fakemass, A,B,C,S, log = True)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S, log = True)
  
    else:
  
      histos = {"dsTau":     selec_C_DsTau[var],
               "dsStarTau": selec_C_DsStarTau[var], 
               "dsMu":     selec_C_DsMu[var], 
               "dsStarMu":  selec_C_DsStarMu[var], 
               "hb"       : selec_C_Hb[var], 
               "data"     : selec_C_Data[var]}
      stackedPlot(histos, var, hb_scale, fakemass,A,B,C,S, mlow, mhigh)
      #stackedPlot(histos, var, hb_scale, fakemass,A,B,C,S, log = True)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S, log = True)

    """ 
    #also do it for weighted_reco
    if 'reco_1' in var: #recos are available
      name = var[0:-7] 
      histos = {"dsMu":     getWeightedReco(selec_S_DsMu,      name, selec, norm = False), 
               "dsTau":     getWeightedReco(selec_S_DsTau,     name, selec, norm = False),
               "dsStarMu":  getWeightedReco(selec_S_DsStarMu,  name, selec, norm = False),
               "dsStarTau": getWeightedReco(selec_S_DsStarTau, name, selec, norm = False),
               "hb"       : selec_S_Hb[var], 
               "combL"    : selecDataL[var], 
               "combR"    : selecDataR[var], 
               "data"     : selec_S_Data[var]}
      stackedPlot(histos, name + "_weighted_reco", hb_scale, fakemass, A,B,C,S, mlow, mhigh)
      #stackedPlot(histos, name + "_weighted_reco", hb_scale, fakemass, A,B,C,S, log = True)
      #normPlot(dict(list(histos.items())[:-1]), name + "_weighted_reco", fakemass, A,B,C,S)
      #normPlot(dict(list(histos.items())[:-1]), name + "_weighted_reco", fakemass, A,B,C,S, log = True)
    """ 
  """

  #include Gen!
  selec_S_DsMu      = createHistos(selecDsMu      + signalRegion, rdfSig)

  #Method comparison
  prefix = ["e_star", "q2","m2_miss","cosMuW","cosPhiDs","cosPlaneBs","bs_pt", "e_star"]
  for name in prefix:
    print(f" ===> Compare Bs reco methods for {name}")
    # for all prefixes get the variables, f.e. m2_miss_gen, m2_miss_coll , ...
    histos = {var: selec_S_DsMu[var] for var in models.keys() if name in var}
    methodDistinction(histos,name,selecDsMu      + signalRegion)
  """

createPlots(baseline)
#createPlots(bkg)
#createPlots(baseline + '&' + addOn1 + '&' + addOn2)

