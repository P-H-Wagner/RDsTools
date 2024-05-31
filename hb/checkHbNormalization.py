import ROOT
import argparse
import os

import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/plots"))

from helper import bsMass_
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

sig   = "22_05_2024_15_50_34"                      #"26_04_2024_16_28_22" #old MC 
hb    = "22_05_2024_16_59_43" #"10_04_2024_00_32_53"#"09_04_2024_18_44_36"  #"22_03_2024_17_34_53" #old inclusive
data  = "26_04_2024_18_08_15" #data -> apply SB method
bs    = "22_05_2024_13_08_32"
bplus = "22_05_2024_13_08_21"
b0    = "22_05_2024_13_08_08"
lambdab    = "22_05_2024_23_22_26"

nSignalRegion    = 3 #how many sigma until the signal region stops
nSidebands     = 5 #how many sigma until the sb starts
sbWidth = 1 # sb width

dsMass_ = 1.96834
bsMass_ = 5.36688
phiMass_ = 1.019461
#########################################
## SELECTIONS                          ##
#########################################

## BASIC

selBasic        = f"((dsMu_m < {bsMass_}) && (k1_charge*k2_charge < 0) && (mu_charge*pi_charge <0) && (gen_match_success == 1))"
selBasicHb      = selBasic + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) "    # exlude signals from hb

#selDsMu         = selBasic + " && (gen_sig == 0)"                                              # select Ds  Mu 
#selDsTau        = selBasic + " && (gen_sig == 1)"                                              # select Ds  Tau 
#selDsStarMu     = selBasic + " && (gen_sig == 10)"                                              # select Ds* Mu
#selDsStarTau    = selBasic + " && (gen_sig == 11)"                                              # select Ds* Tau
#
#selMu           = selBasic + " && (gen_sig == 0 || gen_sig == 10)"                                  # select mu  signals
#selTau          = selBasic + " && (gen_sig == 1 || gen_sig == 11)"                                  # select tau signals
#selDs           = selBasic + " && (gen_sig == 0 || gen_sig == 1)"                                  # select Ds signals
#selDsStar       = selBasic + " && (gen_sig == 10 || gen_sig == 11)"                                  # select Ds star signals

def getRdf(dateTime, debug = False):

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
      #color,_ = getColorAndLabel(var)
      #histos[var].SetLineColor(color)
      histos[var].SetLineWidth(linewidth)
  print(f"      Selection: {selection} DONE")
  return histos

def getColorAndLabelSignalDistinction(key):

  "Adapt colors for multi histogram plot, f.e. mu vs tau"

  labels = []
  colors = []

  if"stacked" in key:
    colors = {}
   
    labels.append(r"Bs"); labels.append(r"B0"); labels.append(r"B+");   labels.append(r"Hb"); labels.append(r"Lambda");
    colors = {"b0":   ROOT.kRed,
              "bplus":  ROOT.kGreen,
              "bs":   ROOT.kBlue,
              "lambdab":         ROOT.kOrange,
              "hb":         ROOT.kBlack}

  if "hb" in key:
    hbColor,hbLabel = getHbColorAndLabel()
    labels.append(hbLabel)
    colors.append(hbColor)

  return (colors,labels)


def stackedPlot(histos, var, scaleB0, scaleBplus, scaleLambdab, log = False):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    print("adapting fill style..", key)
    histos[key].SetFillColor(color[key])
    histos[key].SetLineColor(color[key])
    if key == "hb":
       histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["b0"].GetNbinsX()
  start = histos["b0"].GetXaxis().GetBinLowEdge(1)   # first bin
  stop  = histos["b0"].GetXaxis().GetBinUpEdge(bins) # last bin


  # stacked histo
  hs = ROOT.THStack(var,"")
  # error histo
  hErr = ROOT.TH1D(var,"",bins,start,stop)
  hErr.SetName("err")
  hErr.GetXaxis().SetTitleOffset(1.3)
  hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())  
  hErr.SetLineWidth(0)
  # special error fillstyle 
  hErr.SetFillColor(ROOT.kGray+1)
  hErr.SetFillStyle(3344)

  #scale hb to other mc

  histos["b0"].Scale(scaleB0)
  histos["bplus"].Scale(scaleBplus)
  histos["lambdab"].Scale(scaleLambdab)

  nHb = histos["b0"].Integral() + histos["bplus"].Integral() + histos["bs"].Integral() + histos["lambdab"].Integral() 

  ###################################
  ## Scale MC to data              ##
  ###################################
  for key in histos.keys():
    if ('hb' not in key) :

      histos[key].Scale((histos["hb"].Integral() / nHb))
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
  leg.SetNColumns(3)
  leg.AddEntry(histos["bplus"]      ,'B_{+}#rightarrow Double C' ,'F' )
  leg.AddEntry(histos["bs"]         ,'B_{s}#rightarrow Double C'  ,'F' )
  leg.AddEntry(histos["b0"]         ,'B_{0}#rightarrow Double C' ,'F' )
  leg.AddEntry(histos["hb"]         ,'Inclusive Hb','LEP')
  leg.AddEntry(histos["lambdab"]  ,'#Lambda_{b}#rightarrow Double C','F' )
  leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
  #plot mainpad
  if log:
    ROOT.gPad.SetLogy()
    histos["hb"].GetYaxis().SetTitle('Events')
    histos["hb"].GetYaxis().SetRangeUser(1e-3, hs.GetMaximum()*1000)
    histos["hb"].SetMinimum(1.001) # avoid idsplaying tons of comb bkg
  else:
    histos["hb"].GetYaxis().SetRangeUser(1e-3, histos["hb"].GetMaximum()*1.5)
    histos["hb"].GetYaxis().SetTitle('Events')
  
  #draw data with uncertainty
  histos["hb"].Draw("EP")
  #draw stacked histo
  hs.Draw("HIST SAME")
  #draw data again (used for plotting uses)
  histos["hb"].Draw("EP SAME")
  #draw stat uncertainty of MC, (can not be directly done with a stacked histo)
  hErr.Draw("E2 SAME")
  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = 'X fb^{-1}')
  
  #plot ratiopad
  ratio_pad.cd()
  ratio = histos["hb"].Clone()
  ratio.SetName(ratio.GetName()+'_ratio')
  ratio.Divide(hErr)
  ratio_stats = histos["hb"].Clone()
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
    if ('hb' not in key) :
      h = histos[key].Clone()
      h.Divide(hErr)
      norm_stack.Add(h)
  
  norm_stack.Draw('hist same')

  #print(f"I count {hComb.Integral()} comb. bkg. events for the variable {var}")
 
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

  toSave = "/work/pahwagne/RDsTools/hb/checkNorm/"
  if not os.path.exists(toSave):
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")

  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/hb/checkNorm/log/{var}.pdf")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/hb/checkNorm/{var}.pdf")
  print(f"===> Produced plot: {var}.pdf")



# Create rdf from tree
print(f" ===> Start creating RDataFrames")
chainSig,     rdfSig     = getRdf(sig)
chainHb,      rdfHb      = getRdf(hb)#,   debug = True)
chainB0,      rdfB0      = getRdf(b0)#,   debug = True)
chainBplus,   rdfBplus   = getRdf(bplus)#,   debug = True)
chainBs,      rdfBs      = getRdf(bs)#,   debug = True)
chainLambdab, rdfLambdab = getRdf(lambdab)#,   debug = True)
#chainData,rdfData        = getRdf(data)#skimmed = "baseline") #pick already skimmed file!

# get scale for each component
print(selBasicHb)

lumiBs     = chainBs.GetEntries(selBasicHb + " && (gen_sig == 315)") / 0.0144 

lumiB0     = chainB0.GetEntries(selBasicHb + " && (gen_sig == 213)") / 0.0177
scaleB0    = lumiBs/lumiB0

lumiBplus  = chainBplus.GetEntries(selBasicHb + " && (gen_sig == 114)") / 0.0171
scaleBplus = lumiBs/lumiBplus

lumiLambdab  = chainLambdab.GetEntries(selBasicHb + " && (gen_sig == 416)") / 0.04400
scaleLambdab = lumiBs/lumiLambdab

#for all variables except mass plot only signal region (indicated with 'S')
selec_Hb        = createHistos(selBasicHb          , rdfHb , gen = False)
selec_B0        = createHistos(selBasicHb          , rdfB0 , gen = False)
selec_Bs        = createHistos(selBasicHb          , rdfBs , gen = False)
selec_Bplus     = createHistos(selBasicHb       , rdfBplus , gen = False)
selec_Lambdab   = createHistos(selBasicHb     , rdfLambdab , gen = False)


for var in models.keys():
  print(f"===> Producing stacked plot for variable: {var}") 
  if "gen" in var:
    #skip gen variables
    print("This is a gen variable... skip!")
    continue

  if var != "phiPi_m":

    histos = {"b0":     selec_B0[var],
             "bs":      selec_Bs[var], 
             "bplus":   selec_Bplus[var], 
             "lambdab": selec_Lambdab[var], 
             "hb"     : selec_Hb[var]} 
    stackedPlot(histos, var, scaleB0, scaleBplus, scaleLambdab)
