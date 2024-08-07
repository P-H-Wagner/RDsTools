import ROOT
import argparse
import os
import sys

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from sidebands import getSigma, getABCS
from helper import * 
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
ROOT.TH1.SetDefaultSumw2() #apply weights!

#########################################
## INPUT                               ##
#########################################


#########################################
## SELECTIONS                          ##
#########################################

## BASIC

bkg_enhancing_1    = f"(mu_rel_iso_03 > 0.3) && (lxy_ds_sig < 10) && (dsMu_m < 5.366)" #models the comb
bkg_enhancing_2    = f"(dsMu_m > 5.366)" #models the comb

sign_flip_mu_pi    = f"(pi_charge * mu_charge > 0)"
sign_flip_kk       = f"(k1_charge * k2_charge > 0)"
sign_flip_both     = f"((pi_charge * mu_charge > 0) || (k1_charge * k2_charge > 0))"


#########################################
## CREATE RDF FROM TREE                ##
#########################################

def getRdf(dateTimes, debug = None, skimmed = None, pastNN = None):
 
 
  chain = ROOT.TChain("tree")

  if ((not pastNN) and (not isinstance(dateTimes, list))):
    print("dateTimes must be a list of strings")

  if pastNN:

    print(f"picking past NN file ...") #only one per channel, already chained!
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{dateTimes}.root"
    n = chain.Add(files)
    if n > 1: print("alert, finding more than one past NN tree for this channel!")
    import pdb
    rdf = ROOT.RDataFrame(chain)
    return (chain,rdf)


  if debug:
    print(f"Picking {debug} file(s) for debugging ...")
    fileList = os.listdir(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/")
    #files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/" +  fileList[0] # pick just the first one 
    for i in range(debug):
      try: chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/" +  fileList[i])
      except: chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTimes[0]}/*")

    rdf = ROOT.RDataFrame(chain)
    return (chain,rdf)


  for dateTime in dateTimes:

    if skimmed:
      files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_{skimmed}_{dateTime}.root" #data skimmed with selection 'skimmed'
      #files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_bkg_{dateTime}.root"  # data skimmed with kkpimu > Bs for closure
      print(f"Appending {files}")

    else:
      #access the flat ntuples
      files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test

    #chain them all
    chain.Add(files)

  #create rdf from tree
  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)
 


# Create rdf from tree
print(f" ===> Start creating RDataFrames")
chainData_unc,   rdfData_unc        = getRdf(data_unc[0:1],  debug = False)  
chainData_cons,  rdfData_cons       = getRdf(data_cons[0:1], debug = False)  
chainSigSB,      rdfSigSB           = getRdf(sig_unc,   debug = False) # FOR SB FIT


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


def getReco1Scale(name, selection, start, stop, chain = chainSigSB):

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
  chainSigSB.Project("muNew",    "phiPi_m", selMu + f" & (phiPi_m > {start} ) & (phiPi_m < {stop} )" )

  if (pureHb.Integral() == 0): 
    #avoid dividing by zero (probably muons are also empty!)
    scale_hb = 0
  else:
    scale_hb = muNew.Integral() * pureHb.Integral() / pureMu.Integral()
  print( "#Hb events for MC sample is:", scale_hb)
  return scale_hb

def normPlot(histos, var, fakemass, A,B,C,S, region, log = False, method = ""):

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(ROOT.kRed)
    histos[key].SetLineWidth(3)

  # take as model
  #bins  = histos["dsMu"].GetNbinsX()
  #start = histos["dsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  #stop  = histos["dsMu"].GetXaxis().GetBinUpEdge(bins) # last bin

  # mass histo for fakemass
  #hComb = ROOT.TH1D("mass","mass",bins,start,stop)
  #hComb.GetXaxis().SetTitleOffset(1.3)
  #hComb.GetYaxis().SetTitleSize(1*hComb.GetYaxis().GetTitleSize())
  #hComb.GetXaxis().SetTitleSize(1*hComb.GetXaxis().GetTitleSize())  
  #hComb.SetName("mass")
  #hComb.SetFillStyle(0)
  #hComb.SetLineColor(ROOT.kGray + 2)
  #hComb.SetLineWidth(3)

  ###################################
  ## Combinatorial treatment       ##
  ###################################

  if var != "phiPi_m" and "combL" in histos.keys():
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
    hComb = histos["data_sf"].Clone() #sign flip is modeling comb

  hComb.SetFillStyle(0)
  hComb.SetLineColor(ROOT.kBlue)
  hComb.SetLineWidth(3)

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

  if "combL" in histos.keys():
    label = 'Sideband method'

  else:
    label = 'Sign Flip'
 
  #legend
  leg = ROOT.TLegend(.2,.72,.88,0.88)
  leg.SetBorderSize(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetTextFont(42)
  leg.SetTextSize(0.03)
  leg.SetNColumns(1)
  leg.AddEntry(histos["data"]   ,'Data with bkg. enhancing sel.' ,'L' )
  leg.AddEntry(hComb                ,label  ,'L' )
  
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

  if "combL" in histos.keys():
    name = 'sideband'

  else:
    name = 'signflip'


  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/comb/closure/log/{var}_{name}_{region}_{method}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/comb/closure/log/{var}_{name}_{region}_{method}.png")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/comb/closure/{var}_{name}_{region}_{method}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/comb/closure/{var}_{name}_{region}_{method}.png")
  print(f"===> Produced plot: {var}.pdf")


def createPlots(selec, sign_flip):

  if selec == bkg_enhancing_1: region = "low_mass"
  if selec == bkg_enhancing_2: region = "high_mass"

  #for ds mu mass peak fix selection! To discuss with riccardo.. which one?
  # for now lets pick baseline...
  sigma, h          = getSigma(rdfSigSB, "phiPi_m", baseline + "& gen_sig == 0")

  print("here")
  A, B, C, S        = getABCS( rdfData_unc, selec , "phiPi_m", sigma, h, binsFake = 21, nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)
  #A = 24158.821588264094; B = 226543.39752377832; C = 21219.183790935822; S = 19829.186060030886; #use for debugging (faster) 
  print("done") 
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

  # unconstrained data, model using sideband method
  selec_S_Data_unc      = createHistos(selec          + signalRegion, rdfData_unc , gen = False)
  selec_S_DataL_unc     = createHistos(selec          + leftSB, rdfData_unc ,       gen = False)
  selec_S_DataR_unc     = createHistos(selec          + rightSB, rdfData_unc ,      gen = False)

  #constrained data, model with sign flip
  selec_S_Data_c_bh     = createHistos(selec          + signalRegion, rdfData_cons ,   gen = False)
  selec_S_Data_c_sf     = createHistos(sign_flip       + signalRegion, rdfData_cons ,   gen = False)
 

  #get label for signflip method
  if ("mu" in sign_flip)     and ("k1" not in sign_flip): method = "mu_pi" 
  if ("mu" not in sign_flip) and ("k1" in sign_flip): method  = "kk" 
  else: method = "both"
 
  for var in models.keys():

    print(f"===> Producing stacked plot for variable: {var}") 
    if "gen" in var:
      #skip gen variables
      print("This is a gen variable... skip!")
      continue
  
    if var != "phiPi_m":
  
      histos = {
               "combL"    : selec_S_DataL_unc[var], 
               "combR"    : selec_S_DataR_unc[var], 
               "data"     : selec_S_Data_unc[var]}
      
      
      normPlot(dict(list(histos.items())), var, fakemass, A,B,C,S, region = region, method = method)


      histos = {
               "data"     : selec_S_Data_c_bh[var], 
               "data_sf"     : selec_S_Data_c_sf[var]}


      normPlot(dict(list(histos.items())), var, fakemass, A,B,C,S, region = region, method = method)
 
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

createPlots(bkg_enhancing_1, sign_flip_mu_pi)
createPlots(bkg_enhancing_1, sign_flip_kk)
createPlots(bkg_enhancing_1, sign_flip_both)
#createPlots(bkg_enhancing_2)

