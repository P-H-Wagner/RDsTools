import ROOT
import argparse
import os
import sys

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio
from helper import * 
from histModels import models, pastNN_models, pastNN_2Dmodels

import numpy as np
from cms_style import CMS_lumi

def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'

# parsing
parser = argparse.ArgumentParser()
parser.add_argument('constrained', type=boolean_string) #use constrained sampels?
parser.add_argument('newHb',  type=boolean_string) #use new hb?
parser.add_argument('pastNN', type=boolean_string) #use samples after NN?
args = parser.parse_args()

print(f"====> Running constrained fit? {args.constrained}")
if args.pastNN: models.update(pastNN_models)


# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

##############################
# Load chain into RDataFrame #
##############################



def getRdf(dateTimes, debug = None, skimmed = None, pastNN = None):

  print(dateTimes)
  chain = ROOT.TChain("tree")

  if ((not pastNN) and (not isinstance(dateTimes, list))):
    print("dateTimes must be a list of strings")

  if pastNN:
    print(f"picking past NN files ...")
    #dateTimes here is a string like: "data_26Sep..._cons"
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{dateTimes}/*"
    print(files)
    n = chain.Add(files)
    print(f"Chaining {n} files for this channel")

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


########################################
# Assign selections depending on the   #
# MC signal ID                         #
########################################

class selections:

  def __init__(self, selec):
    self.hb =        selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1 ) "
    self.bs =        selec + " && (300 <= gen_sig) && (gen_sig < 400) "
    self.b0 =        selec + " && (200 <= gen_sig) && (gen_sig < 300) "
    self.bplus =     selec + " && (100 <= gen_sig) && (gen_sig < 200) "
    self.dsMu =      selec + " && (gen_sig == 0)" 
    self.dsTau =     selec + " && (gen_sig == 1)"
    self.dsStarMu =  selec + " && (gen_sig == 10)" 
    self.dsStarTau = selec + " && (gen_sig == 11)"
    self.bare = selec

#######################################
# Define signal regions and sidebands #
#######################################

def getRegions(sigma):

  #signal region
  mlow   = dsMass_ - nSignalRegion*sigma
  mhigh  = dsMass_ + nSignalRegion*sigma

  #sideband start
  mlow2  = dsMass_ - nSidebands*sigma
  mhigh2 = dsMass_ + nSidebands*sigma

  #sideband stops
  mlow3  = mlow2  - sbWidth*sigma
  mhigh3 = mhigh2 + sbWidth*sigma
  
  signalRegion      = f"&& ({mlow} < phiPi_m) && (phiPi_m < {mhigh})"
  anti_signalRegion = f"&& ((({mlow3} < phiPi_m) && (phiPi_m < {mlow2})) || (({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3}))) "
  leftSB            = f"&& ({mlow3} < phiPi_m) && (phiPi_m < {mlow2})"
  rightSB           = f"&& ({mhigh2} < phiPi_m) && (phiPi_m < {mhigh3})"

  return mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, leftSB, rightSB

#########################################
## INPUT SAMPLES                       ##
#########################################

# Create rdf from tree
print(f" ===> Start creating RDataFrames")


if args.pastNN:

  print("Producing post-NN plots")

  #save the NN name
  import re
  pattern = r"\d{2}[A-Za-z]{3}\d{4}_\d{2}h\d{2}m\d{2}s"
  match = re.search(pattern,sig_cons_pastNN)
  nnModel = match.group(0)



  if args.constrained:
    print("constrained data after NN not processed yet..")
    chainSigSB, rdfSigSB     = getRdf(sig_cons_pastNN       , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   ) # FOR SB FIT
    chainSig,   rdfSig       = getRdf(sig_cons_pastNN       , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   )#,  debug = True)
    chainHb,    rdfHb        = getRdf(hb_cons_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBs,    rdfBs        = getRdf(bs_cons_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainB0,    rdfB0        = getRdf(b0_cons_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBplus, rdfBplus     = getRdf(bplus_cons_pastNN     , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainData,  rdfData      = getRdf(data_cons_pastNN      , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv")  #skimmed = "baseline") #pick already skimmed file!



  else:
    chainSigSB, rdfSigSB     = getRdf(sig_unc_pastNN       , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   ) # FOR SB FIT
    chainSig,   rdfSig       = getRdf(sig_unc_pastNN       , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   )#,  debug = True)
    chainHb,    rdfHb        = getRdf(hb_unc_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBs,    rdfBs        = getRdf(bs_unc_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainB0,    rdfB0        = getRdf(b0_unc_pastNN        , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBplus, rdfBplus     = getRdf(bplus_unc_pastNN     , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainData,  rdfData      = getRdf(data_unc_pastNN      , pastNN = args.pastNN )#, debug = 10)#,     skimmed = "ma_cut_wout_fv")  #skimmed = "baseline") #pick already skimmed file!


else:

  print("Producing pre-NN plots")
  if args.constrained:
    chainSigSB, rdfSigSB     = getRdf(sig_cons      , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"   ) # FOR SB FIT
    chainSig, rdfSig         = getRdf(sig_cons      , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"   )#,  debug = True)
    chainHb,  rdfHb          = getRdf(hb_cons       , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBs,  rdfBs          = getRdf(bs_cons       , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainB0,  rdfB0          = getRdf(b0_cons       , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBplus,  rdfBplus    = getRdf(bplus_cons    , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainData,rdfData        = getRdf(data_cons     , skimmed = "base_wout_tv")#, debug = 1)#,      skimmed = "ma_cut_wout_fv")  #skimmed = "baseline") #pick already skimmed file!
  else:
  
    chainSigSB, rdfSigSB     = getRdf(sig_unc       , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   ) # FOR SB FIT
    chainSig,   rdfSig       = getRdf(sig_unc       , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"   )#,  debug = True)
    chainHb,    rdfHb        = getRdf(hb_unc        , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBs,    rdfBs        = getRdf(bs_unc        , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainB0,    rdfB0        = getRdf(b0_unc        , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainBplus, rdfBplus     = getRdf(bplus_unc     , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv"    )#,   debug = True)
    chainData,  rdfData      = getRdf(data_unc      , skimmed = "base_wout_tv" )#, debug = 10)#,     skimmed = "ma_cut_wout_fv")  #skimmed = "baseline") #pick already skimmed file!
  
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
              "bs":         ROOT.kRed - 9,
              "b0":         ROOT.kRed - 6,
              "bplus":      ROOT.kRed - 5,
              "combL":      ROOT.kGray+1,
              "combR":      ROOT.kGray+1,
              "data_sf":    ROOT.kGray+1,
              "data_sf_kk":    ROOT.kGray+1,
              "data_sf_pimu":    ROOT.kGray+1,
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
    if norm: 
      try: hist.Scale(1/hist.Integral())
      except: pass 
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

def createHistos(selection,rdf, linewidth = 2, gen = True, massPlot = False, variables = None):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"
  "If <variables> are given, only these variables from <models> are created" #f.e. for the phiPi mass or DsMu with different selection

  histos = {}
  for var,model in models.items(): 
      if variables:
        if var not in variables: 
          print(f"cannot find {var}");     
          continue 

      if "gen" in var and not gen:
        #skip gen variables
        print("This is a gen variable... skip!")
        continue

      tofill = var
      print("i am here")
      if (massPlot == True and var != "phiPi_m"):
        continue

      elif (massPlot == True and var == "phiPi_m"):
        tofill = f"(run==1) * 0.9995* phiPi_m + (run!=1) * phiPi_m"
        print("correcting and plotting phiPi mass..")
        histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], "m_corr")

      else:
        print("Filling for variable ", var)
        histos[var] = rdf.Filter(selection).Histo1D(model[0], tofill)

      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetMaximum())
      
      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(linewidth)


  print(f"      Selection: {selection} DONE")

  return histos

def create2DHistos(selection,rdf, linewidth = 2, gen = True, massPlot = False):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"

  histos2D = {}
   
  for var,model in pastNN_2Dmodels.items(): 

      var1 = var.split("_")[0]
      var2 = var.split("_")[1]

      print("Filling for variable ", var1, " x " , var2)
      histos2D[var] = rdf.Filter(selection).Histo2D(model[0], var1, var2)

      histos2D[var].GetXaxis().SetTitle(model[1])
      histos2D[var].SetMaximum(1.2 * histos2D[var].GetMaximum())
      
      #get default colors 
      color,_ = getColorAndLabel(var)
      histos2D[var].SetLineColor(color)
      histos2D[var].SetLineWidth(linewidth)


  print(f"      Selection: {selection} DONE")

  return histos2D

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

  return scale_hb


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



def stackedPlot(histos, var, hb_scale, scale_kk, scale_bkg, scale_n, mlow, mhigh, constrained, log = False, fakemass = None, A = None, B = None, C = None, S = None, bs = None, b0 = None, bplus = None, newHb = False):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

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
  hErr.SetFillColor(ROOT.kBlack)
  hErr.SetFillStyle(3244)

  #scale hb to other mc

  if hb_scale == 0: 
    histos["hb"].Scale(0.0 )
  else:
    histos["hb"].Scale(hb_scale / histos["hb"].Integral())

  if newHb:
    histos["bs"].Scale(histos["hb"].Integral()    * bs/ histos["bs"].Integral())
    histos["b0"].Scale(histos["hb"].Integral()    * b0/ histos["b0"].Integral())
    histos["bplus"].Scale(histos["hb"].Integral() * bplus/ histos["bplus"].Integral())
    #change from inclusive to exlusive
    histos["hb"] = histos["bs"].Clone()
    histos["hb"].Add(histos["b0"])
    histos["hb"].Add(histos["bplus"])


  nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()
  print("==========> number of signals inside stacked", nSig)


  ###################################
  ## Combinatorial treatment       ##
  ###################################

  if constrained:
    #simple: take sign slipped data as comb
    #hComb = histos["data_sf"].Clone()
    #hComb.Scale(0.1)

    hComb      = histos["data_sf_kk"].Clone()
    hComb.Scale(1.0/hComb.Integral())

    hComb_pimu = histos["data_sf_pimu"].Clone()
    hComb_pimu.Scale(1.0/hComb_pimu.Integral())

    hComb.Scale(scale_kk)
    hComb_pimu.Scale(1.0 - scale_kk)
    hComb.Add(hComb_pimu)

    hComb.Scale(scale_bkg/hComb.Integral())
    print("sclae bkg is: ", scale_bkg)
    print("hcomb has now integral: ", hComb.Integral())

  else: 
    #do sideband method
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
      
  nComb = hComb.Integral()
  
  histos["comb"] = hComb # append to create datacard later
  ###################################
  ## Scale MC to data              ##
  ###################################


  if constrained:

    #histos['comb'].Scale( histos["data"].Integral()  )

    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    for key in histos.keys():

      #normalize and scale the signals and hb
      if ('comb' not in key) and ('data' not in key):
        histos[key].Scale(1.0/nSig)
        histos[key].Scale(1.0 - scale_bkg)

    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    print("nsig has integral:", nSig)

    for key in histos.keys():
      #normalize to data!
      if ('data' not in key):
        histos[key].Scale( histos["data"].Integral() / (nSig + nComb))
        #print("scale_n inside :", histos["data"].Integral())
        print("rescaling...", key)
  else:
    #For sideband method, combinatorial is already normalized to data!
    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key):
        histos[key].Scale((histos["data"].Integral() - nComb) / nSig)

  # add comb to stack and err 
  hErr.Add(hComb)
  hs.Add(hComb)


  if newHb:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key):

        hErr.Add(histos[key])
        hs.Add(histos[key])
  else:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key):
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
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','F' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'F' )

  if newHb:
    leg.AddEntry(histos["bs"]       ,'B_{s}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["b0"]       ,'B_{0}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["bplus"]    ,'B_{+}#rightarrow D_{s} + #mu ' ,'F' )
  else:
    leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  leg.AddEntry(histos["data"]       ,'Data','LEP')
  leg.AddEntry(histos["comb"]                           ,'Comb. Bkg.'  ,'F' )
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
  #hErr.Draw("E2 SAME")

  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  
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

  if newHb:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key):
        h = histos[key].Clone()
        h.Divide(hErr)
        norm_stack.Add(h)

  else: 
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key):
        h = histos[key].Clone()
        h.Divide(hErr)
        norm_stack.Add(h)

  hDummy = histos["comb"].Clone() #do it also for the comb
  hDummy.Divide(hErr)
  norm_stack.Add(hDummy)
  
  norm_stack.Draw('hist same')

  #subplot
  line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1., ratio.GetXaxis().GetXmax(), 1.)
  line.SetLineColor(ROOT.kBlack)
  line.SetLineWidth(1)
  ratio_stats.GetYaxis().CenterTitle()
  
  ratio_stats.Draw('E2')
  norm_stack.Draw('hist same')
  ratio_stats.Draw('E2 same')
  line.Draw('same')
  ratio.Draw('EP same')
  ROOT.gPad.RedrawAxis()
  
  #saving
  c1.Modified()
  c1.Update()



  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if args.pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  toSave = f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}/log/{var}_{name}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}/log/{var}_{name}.png")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}/{var}_{name}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")


  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()

  return returnHisto


def stacked2DPlot(histos, var, hb_scale, mlow, mhigh, constrained, log = False, fakemass = None, A = None, B = None, C = None, S = None, bs = None, b0 = None, bplus = None, newHb = False):

  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillColorAlpha(color[key],0.3)
    histos[key].SetLineColor(color[key])


    if key == "data":
       histos[key].SetMarkerStyle(8) 

  # take as model
  #bins  = histos["dsMu"].GetNbinsX()
  #start = histos["dsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  #stop  = histos["dsMu"].GetXaxis().GetBinUpEdge(bins) # last bin


  # stacked histo
  hs = ROOT.THStack(var,"")
  # error histo
  #hErr = ROOT.TH1D(var,"",bins,start,stop)
  #hErr.SetName("err")
  #hErr.GetXaxis().SetTitleOffset(1.3)
  #hErr.GetYaxis().SetTitleSize(1*hErr.GetYaxis().GetTitleSize())
  #hErr.GetXaxis().SetTitleSize(1*hErr.GetXaxis().GetTitleSize())  
  #hErr.SetLineWidth(0)
  # mass histo for fakemass
  #hComb = hErr.Clone()
  #hComb.SetName("mass")
  #hComb.SetFillColor(ROOT.kGray)
  #hComb.SetLineColor(ROOT.kGray)
  # special error fillstyle 
  #hErr.SetFillColor(ROOT.kGray+1)
  #hErr.SetFillStyle(3344)

  #scale hb to other mc

  print("Scale of Hb: ", hb_scale )
  if hb_scale == 0: 
    histos["hb"].Scale(0.0 )
  else:
    histos["hb"].Scale(hb_scale / histos["hb"].Integral())

  if newHb:
    histos["bs"].Scale(histos["hb"].Integral()    * bs/ histos["bs"].Integral())
    histos["b0"].Scale(histos["hb"].Integral()    * b0/ histos["b0"].Integral())
    histos["bplus"].Scale(histos["hb"].Integral() * bplus/ histos["bplus"].Integral())
    #change from inclusive to exlusive
    histos["hb"] = histos["bs"].Clone()
    histos["hb"].Add(histos["b0"])
    histos["hb"].Add(histos["bplus"])


  nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

  ###################################
  ## Combinatorial treatment       ##
  ###################################

  if constrained:
    hComb = histos["data_sf"].Clone()

  else: 
    #do sideband method
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
      
  nComb = hComb.Integral()
  
  hComb.SetName("mass")
  hComb.SetFillColorAlpha(ROOT.kGray, 0.2)
  hComb.SetLineColor(ROOT.kGray)

  # add comb to stack and err 
  #hErr.Add(hComb)
  hs.Add(hComb)
  
  histos["comb"] = hComb # append to create datacard later
  ###################################
  ## Scale MC to data              ##
  ###################################


  for key in histos.keys():
    if ('comb' not in key) and ('data' not in key):

      histos[key].Scale((histos["data"].Integral() - nComb) / nSig)

  if newHb:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key):

        #hErr.Add(histos[key])
        hs.Add(histos[key])
  else:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key):
        #hErr.Add(histos[key])
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
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','F' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'F' )

  if newHb:
    leg.AddEntry(histos["bs"]       ,'B_{s}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["b0"]       ,'B_{0}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["bplus"]    ,'B_{+}#rightarrow D_{s} + #mu ' ,'F' )
  else:
    leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
  leg.AddEntry(histos["data"]       ,'Data','LEP')
  leg.AddEntry(hComb                           ,'Comb. Bkg.'  ,'F' )
  #leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
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
  histos["data"].Draw("LEGO1")
  #draw stacked histo
  hs.Draw("LEGO1 SAME")
  #draw data again (used for plotting uses)
  histos["data"].Draw("LEGO1 EP SAME")
  #draw stat uncertainty of MC, (can not be directly done with a stacked histo)
  #hErr.Draw("E2 SAME")
  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  #ROOT.gPad.SetPhi(90)
 
  """ 
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

  if newHb:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key):
        h = histos[key].Clone()
        h.Divide(hErr)
        norm_stack.Add(h)

  else: 
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key):
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

  """
  toSave = "/work/pahwagne/RDsTools/plots/cmsplots/2D/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if args.pastNN: name += "_pastNN"
  else: name += ""


  c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/2D/{var}_{name}.pdf")
  c1.SaveAs(f"/work/pahwagne/RDsTools/plots/cmsplots/2D/{var}_{name}.png")
  print(f"===> Produced plot: {var}.pdf")


  # return a copy, otherwise the histosScaled get overwritten when we change histos in f.e. normplots! (Weird???)
  returnHisto = {}
  for key in histos.keys():
    returnHisto[key] = histos[key].Clone()

  return returnHisto

def normPlot(histos, var, scale_kk, constrained, fakemass = None , A = None ,B = None ,C = None ,S = None ,log = False, newHb = False):

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

  if constrained:
    #simple: take sign slipped data as comb
    hComb      = histos["data_sf_kk"].Clone()
    hComb_pimu = histos["data_sf_pimu"].Clone()

    hComb.Scale(scale_kk)
    hComb_pimu.Scale(1.0 - scale_kk)

    hComb.Add(hComb_pimu)

  else:

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
  leg.AddEntry(histos["dsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'L' )
  leg.AddEntry(histos["dsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'L' )
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','L' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'L' )
  if newHb:
    leg.AddEntry(histos["bs"]       ,'B_{s}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["b0"]       ,'B_{0}#rightarrow D_{s} + #mu ' ,'F' )
    leg.AddEntry(histos["bplus"]    ,'B_{+}#rightarrow D_{s} + #mu ' ,'F' )
  else:
    leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'F' )
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

    print(f"drawing variable {key}")

    if ("comb" in key) or ("hb" in key and newHb == True): 
      #comes later, dont draw sidebands
      continue

    if i == 0:
      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")

  hComb.Draw("HIST SAME") 
  leg.Draw("SAME")
  #ROOT.gPad.RedrawAxis()
  
  CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  #saving
  c1.Modified()
  c1.Update()


  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if args.pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  toSave = f"/work/pahwagne/RDsTools/plots/cmsplots/{folder}"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  


  if log:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{folder}/log/{var}_{name}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{folder}/log/{var}_{name}.png")
  else:
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{folder}/{var}_{name}.pdf")
    c1.SaveAs(f"/work/pahwagne/RDsTools/plots/normplots/{folder}/{var}_{name}.png")
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

def writeDatacard(histos, var, digits = 5, is_2D = None):

  if args.constrained: name = "constrained"
  else: name = "unconstrained"

  if args.newHb: name += "_newHb"
  else: name += ""

  if args.pastNN: name += "_pastNN"
  else: name += ""

  if is_2D: name += "_2D"

  temp = open("/work/pahwagne/RDsTools/fit/datacardTemplate.txt", "rt")
  card = open(f"/work/pahwagne/RDsTools/fit/datacards/datacard_{var}_{name}.txt", "wt")

  # width in datacard template is 17 spaces
  spaces = 22

  dataStr       = str(round(histos["data"].Integral(), digits))
  dataStr      += " "*(spaces - len(dataStr))
  
  dsMuStr       = str(round(histos["dsMu"].Integral(), digits))
  dsMuStr      += " "*(spaces - len(dsMuStr))

  dsTauStr      = str(round(histos["dsTau"].Integral(), digits))
  dsTauStr     += " "*(spaces - len(dsTauStr))

  dsStarMuStr   = str(round(histos["dsStarMu"].Integral(), digits))
  dsStarMuStr  += " "*(spaces - len(dsStarMuStr))

  dsStarTauStr  = str(round(histos["dsStarTau"].Integral(), digits))
  dsStarTauStr += " "*(spaces - len(dsStarTauStr))

  hbStr         = str(round(histos["hb"].Integral(), digits))
  hbStr        += " "*(spaces - len(hbStr))

  combStr       = str(round(histos["comb"].Integral(), digits))  
  #combStr      += " "*(spaces - len(combStr))

  rates = dsMuStr + dsTauStr + dsStarMuStr + dsStarTauStr + hbStr + combStr

  for line in temp:
    if "HOOK_RATES"            in line: line = line.replace("HOOK_RATES", rates )
    if "HOOK_DATA_RATE"      in line: line = line.replace("HOOK_DATA_RATE", dataStr )
    if "HOOK_VAR"            in line: line = line.replace("HOOK_VAR", var )
    if "HOOK_NAME"            in line: line = line.replace("HOOK_NAME", name )
    card.write(line)

  temp.close()
  card.close()

def createPlots(baseline, constrained = False, newHb = False):

  ################################
  # Define important selections  #
  ################################

  high_mass         = f" && (dsMu_m > {bsMass_})"
  low_mass          = f" && (dsMu_m < {bsMass_})"
  wrong_sign        = f" && (((k1_charge*k2_charge > 0) || (mu_charge*pi_charge > 0)))"
  right_sign        = f" && (((k1_charge*k2_charge < 0) && (mu_charge*pi_charge < 0)))"
  kk_wrong          = f" && (k1_charge*k2_charge > 0) && (mu_charge*pi_charge < 0)"
  pimu_wrong        = f" && (mu_charge*pi_charge > 0) && (k1_charge*k2_charge < 0)"
  both_wrong        = f" && (mu_charge*pi_charge > 0) && (k1_charge*k2_charge > 0)"
  

  # Define selections which hold all important selections next to baseline  
  selection_H     = baseline + right_sign             #for Ds + Mu mass, includes high mass region!
  selection       = baseline + right_sign + low_mass  #for all other variables

  # after NN, add score cut to 'selec*'
  score_cut = ""

  if args.pastNN:
 
    score_cut = "&& (score5 <= 0.3)" #&& ((score0 > 0.3) || (score1 > 0.3) || (score2 > 0.5) || (score3 > 0.5))"
    selection_H += score_cut  # && (tv_prob > 0.1)"
    selection   += score_cut  # && (tv_prob > 0.1)"

  #Now define a selection for every MC signal and MC background, slicing with their ID 

  selec_H = selections(selection_H)
  selec   = selections(selection)

  ################################################
  # Extract sigma from fit to MC DsMu signal to  #
  # define signal region.                        #
  ################################################

  sigma_H, h_H          = getSigma(rdfSigSB, "phiPi_m", selec_H.dsMu)
  sigma,   h            = getSigma(rdfSigSB, "phiPi_m", selec.dsMu)

  mlow_H, mhigh_H, mlow2_H, mhigh2_H, mlow3_H, mhigh3_H, signalRegion_H, anti_signalRegion_H, leftSB_H, rightSB_H = getRegions(sigma)
  mlow,   mhigh,   mlow2,   mhigh2,   mlow3,   mhigh3,   signalRegion,   anti_signalRegion,   leftSB,   rightSB   = getRegions(sigma)

  ################################################
  # Fit Data DsMu signal to define A,B,C,D for   #
  # sideband method and fakemass for sideband    #
  ################################################

  if not constrained:

    #do sideband method and extract A, B, C ,D and fakemass histo including high mass region
    A_H, B_H, C_H, S_H  = getABCS( rdfData, selec_H.bare , "phiPi_m", sigma_H, h_H, binsFake = 21, nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)
    #get fakemass histo for the phiPi mass
    fakemass_H = np.genfromtxt('mass_bincontent.csv', delimiter=',')

    #repeat for low mass region only
    A, B, C, S          = getABCS( rdfData, selec.bare, "phiPi_m", sigma,   h,   binsFake = 21, nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)
    fakemass = np.genfromtxt('mass_bincontent.csv', delimiter=',')

  #get proportions from inclusive sample
  hb_tot_H      = rdfHb.Filter(selec_H.hb).Count().GetValue()
  bs_in_hb_H    = rdfHb.Filter(selec_H.bs).Count().GetValue()    / hb_tot_H
  b0_in_hb_H    = rdfHb.Filter(selec_H.b0).Count().GetValue()    / hb_tot_H
  bplus_in_hb_H = rdfHb.Filter(selec_H.bplus).Count().GetValue() / hb_tot_H

  hb_tot        = rdfHb.Filter(selec.hb).Count().GetValue()
  bs_in_hb      = rdfHb.Filter(selec.bs).Count().GetValue()    / hb_tot
  b0_in_hb      = rdfHb.Filter(selec.b0).Count().GetValue()    / hb_tot
  bplus_in_hb   = rdfHb.Filter(selec.bplus).Count().GetValue() / hb_tot

  #create histos returns a dictionary !:)
 
  ###################################################
  # Histograms in signal region (for most variables #
  ###################################################
 
  print("===> filling signal region ...")
  hb_scale_S = getHbScale(selec.hb + signalRegion, selec.dsMu + signalRegion)

  #for all variables except mass plot only signal region (indicated with 'S')
  selec_S_DsMu           = createHistos(selec.dsMu        + signalRegion,    rdfSig        , gen = False)
  selec_S_DsTau          = createHistos(selec.dsTau       + signalRegion,    rdfSig        , gen = False)
  selec_S_DsStarMu       = createHistos(selec.dsStarMu    + signalRegion,    rdfSig        , gen = False)
  selec_S_DsStarTau      = createHistos(selec.dsStarTau   + signalRegion,    rdfSig        , gen = False)
  selec_S_Hb             = createHistos(selec.hb          + signalRegion,    rdfHb         , gen = False)
  selec_S_Data           = createHistos(selec.bare        + signalRegion,    rdfData       , gen = False)

  selec_S_DsMu_2D        = create2DHistos(selec.dsMu      + signalRegion,    rdfSig        , gen = False)
  selec_S_DsTau_2D       = create2DHistos(selec.dsTau     + signalRegion,    rdfSig        , gen = False)
  selec_S_DsStarMu_2D    = create2DHistos(selec.dsStarMu  + signalRegion,    rdfSig        , gen = False)
  selec_S_DsStarTau_2D   = create2DHistos(selec.dsStarTau + signalRegion,    rdfSig        , gen = False)
  selec_S_Hb_2D          = create2DHistos(selec.hb        + signalRegion,    rdfHb         , gen = False)
  selec_S_Data_2D        = create2DHistos(selec.bare      + signalRegion,    rdfData       , gen = False)

  if newHb:
    selec_S_Bs           = createHistos(selec.bs          + signalRegion,    rdfBs         , gen = False)
    selec_S_B0           = createHistos(selec.b0          + signalRegion,    rdfB0         , gen = False)
    selec_S_Bplus        = createHistos(selec.bplus       + signalRegion,    rdfBplus      , gen = False)

    selec_S_Bs_2D        = create2DHistos(selec.bs        + signalRegion,    rdfBs         , gen = False)
    selec_S_B0_2D        = create2DHistos(selec.b0        + signalRegion,    rdfB0         , gen = False)
    selec_S_Bplus_2D     = create2DHistos(selec.bplus     + signalRegion,    rdfBplus      , gen = False)

  if not constrained:
    selec_S_DataL        = createHistos(selec.bare        + leftSB,          rdfData       , gen = False)
    selec_S_DataR        = createHistos(selec.bare        + rightSB,         rdfData       , gen = False)

    selec_S_DataL_2D     = create2DHistos(selec.bare      + leftSB,          rdfData       , gen = False)
    selec_S_DataR_2D     = create2DHistos(selec.bare      + rightSB,         rdfData       , gen = False)

  else: 
    #selec_S_Data_sf      = createHistos(baseline + score_cut + wrong_sign + low_mass + signalRegion,  rdfData       , gen = False)
    selec_S_Data_sf_kk   = createHistos(  baseline + score_cut + kk_wrong   + low_mass + signalRegion,  rdfData       , gen = False)
    selec_S_Data_sf_pimu = createHistos(  baseline + score_cut + pimu_wrong + low_mass + signalRegion,  rdfData       , gen = False)
    selec_S_Data_sf_2D   = create2DHistos(baseline + score_cut + wrong_sign + low_mass + signalRegion,  rdfData       , gen = False)

    #get the signflip scale by fitting the ds mass peak
  
    hbDummy = selec_S_Hb["phiPi_m"].Clone()
    if hb_scale_S == 0:
      hbDummy.Scale(0.0 )
    else:
      hbDummy.Scale(hb_scale_S / hbDummy.Integral())

    hRest    = selec_S_DsMu["phiPi_m"].Clone()
    hRest.Add( selec_S_DsStarMu["phiPi_m"].Clone())
    hRest.Add( selec_S_DsTau["phiPi_m"].Clone())
    hRest.Add( selec_S_DsStarTau["phiPi_m"].Clone())
    hRest.Add( hbDummy)

    global scale_S_kk;
    global scale_S_bkg;
    global scale_S_n;

    print("===> signflip fit...")
    print("We have ",selec_S_Data_sf_kk["phiPi_m"].GetEntries(), "wrong kk sign events")
    print("We have ",selec_S_Data_sf_pimu["phiPi_m"].GetEntries(), "wrong pimu sign events")

    scale_S_kk, scale_S_bkg, scale_S_n = getSignflipRatio(selec_S_Data_sf_kk["phiPi_m"].Clone(),selec_S_Data_sf_pimu["phiPi_m"].Clone(),hRest,selec_S_Data["phiPi_m"].Clone())

  print("===> signal region done...")

  print("===> filling total region ...")
  ###################################################
  # Histograms in complete region (for Ds mass)     #
  ###################################################

  hb_scale_C = getHbScale(selec.hb , selec.dsMu )
  
  # for the ds mass plot we also want to plot the sidebands! (indicated with 'C' for complete)
  selec_C_DsMu           = createHistos(selec.dsMu       ,                  rdfSig        , gen = False, massPlot = True, variables = ["phiPi_m"])
  selec_C_DsTau          = createHistos(selec.dsTau      ,                  rdfSig        , gen = False, massPlot = True, variables = ["phiPi_m"])
  selec_C_DsStarMu       = createHistos(selec.dsStarMu   ,                  rdfSig        , gen = False, massPlot = True, variables = ["phiPi_m"])
  selec_C_DsStarTau      = createHistos(selec.dsStarTau  ,                  rdfSig        , gen = False, massPlot = True, variables = ["phiPi_m"])
  selec_C_Hb             = createHistos(selec.hb         ,                  rdfHb         , gen = False, massPlot = True, variables = ["phiPi_m"])
  selec_C_Data           = createHistos(selec.bare       ,                  rdfData       , gen = False, massPlot = True, variables = ["phiPi_m"])

  selec_C_DsMu_2D        = create2DHistos(selec.dsMu      ,                 rdfSig        , gen = False, massPlot = True )
  selec_C_DsTau_2D       = create2DHistos(selec.dsTau     ,                 rdfSig        , gen = False, massPlot = True )
  selec_C_DsStarMu_2D    = create2DHistos(selec.dsStarMu  ,                 rdfSig        , gen = False, massPlot = True )
  selec_C_DsStarTau_2D   = create2DHistos(selec.dsStarTau ,                 rdfSig        , gen = False, massPlot = True )
  selec_C_Hb_2D          = create2DHistos(selec.hb        ,                 rdfHb         , gen = False, massPlot = True )
  selec_C_Data_2D        = create2DHistos(selec.bare      ,                 rdfData       , gen = False, massPlot = True )


  if newHb:
    selec_C_Bs           = createHistos(selec.bs          ,                 rdfBs         , gen = False, massPlot = True, variables = ["phiPi_m"])
    selec_C_B0           = createHistos(selec.b0          ,                 rdfB0         , gen = False, massPlot = True, variables = ["phiPi_m"])
    selec_C_Bplus        = createHistos(selec.bplus       ,                 rdfBplus      , gen = False, massPlot = True, variables = ["phiPi_m"])

    selec_C_Bs_2D        = create2DHistos(selec.bs        ,                 rdfBs         , gen = False, massPlot = True )
    selec_C_B0_2D        = create2DHistos(selec.b0        ,                 rdfB0         , gen = False, massPlot = True )
    selec_C_Bplus_2D     = create2DHistos(selec.bplus     ,                 rdfBplus      , gen = False, massPlot = True )

  if constrained:
    #selec_C_Data_sf     = createHistos(baseline + score_cut + wrong_sign + low_mass,   rdfData        , gen = False, massPlot = True, variables = ["phiPi_m"])
    selec_C_Data_sf_kk   = createHistos(  baseline + score_cut + kk_wrong   + low_mass,      rdfData       , gen = False, massPlot = True)#, variables = ["phiPi_m","q2_coll"])
    selec_C_Data_sf_pimu = createHistos(  baseline + score_cut + pimu_wrong + low_mass,      rdfData       , gen = False, massPlot = True)#, variables = ["phiPi_m","q2_coll"])
    selec_C_Data_sf_both = createHistos(  baseline + score_cut + both_wrong + low_mass,      rdfData       , gen = False, massPlot = True)#, variables = ["phiPi_m","q2_coll"])
    selec_C_Data_sf_2D   = create2DHistos(baseline + score_cut + wrong_sign + low_mass,      rdfData       , gen = False, massPlot = True)
 

    #get the signflip scale by fitting the ds mass peak
  
    hbDummy = selec_C_Hb["phiPi_m"].Clone()
    if hb_scale_C == 0:
      hbDummy.Scale(0.0 )
    else:
      hbDummy.Scale(hb_scale_C / hbDummy.Integral())

    hRest    = selec_C_DsMu["phiPi_m"].Clone()
    hRest.Add( selec_C_DsStarMu["phiPi_m"].Clone())
    hRest.Add( selec_C_DsTau["phiPi_m"].Clone())
    hRest.Add( selec_C_DsStarTau["phiPi_m"].Clone())
    hRest.Add( hbDummy)

    global scale_C_kk;
    global scale_C_bkg;
    global scale_C_n;


    for var in selec_C_Data_sf_both.keys():

      plotSignFlipShapes(selec_C_Data_sf_kk,selec_C_Data_sf_pimu,selec_C_Data_sf_both,var)

    scale_C_kk, scale_C_bkg, scale_C_n = getSignflipRatio(selec_C_Data_sf_kk["phiPi_m"].Clone(),selec_C_Data_sf_pimu["phiPi_m"].Clone(),hRest,selec_C_Data["phiPi_m"].Clone())
    print("I AM HERE!!")
 
  print("===> total region done...")


  print("===> filling high mass region ...")

  ########################################################
  # Histograms in high&low mass region (For Ds+mu mass)  #
  ########################################################

  hb_scale_H = getHbScale(selec_H.hb , selec_H.dsMu )
  
  # for the ds mass plot we also want to plot the sidebands! (indicated with 'C' for complete)
  selec_H_DsMu           = createHistos(selec_H.dsMu       ,                  rdfSig        , gen = False,  variables = ["phiPi_m", "dsMu_m"])
  selec_H_DsTau          = createHistos(selec_H.dsTau      ,                  rdfSig        , gen = False,  variables = ["phiPi_m", "dsMu_m"])
  selec_H_DsStarMu       = createHistos(selec_H.dsStarMu   ,                  rdfSig        , gen = False,  variables = ["phiPi_m", "dsMu_m"])
  selec_H_DsStarTau      = createHistos(selec_H.dsStarTau  ,                  rdfSig        , gen = False,  variables = ["phiPi_m", "dsMu_m"])
  selec_H_Hb             = createHistos(selec_H.hb         ,                  rdfHb         , gen = False,  variables = ["phiPi_m", "dsMu_m"])
  selec_H_Data           = createHistos(selec_H.bare       ,                  rdfData       , gen = False,  variables = ["phiPi_m", "dsMu_m"])

  selec_H_DsMu_2D        = create2DHistos(selec_H.dsMu      ,                 rdfSig        , gen = False )
  selec_H_DsTau_2D       = create2DHistos(selec_H.dsTau     ,                 rdfSig        , gen = False )
  selec_H_DsStarMu_2D    = create2DHistos(selec_H.dsStarMu  ,                 rdfSig        , gen = False )
  selec_H_DsStarTau_2D   = create2DHistos(selec_H.dsStarTau ,                 rdfSig        , gen = False )
  selec_H_Hb_2D          = create2DHistos(selec_H.hb        ,                 rdfHb         , gen = False )
  selec_H_Data_2D        = create2DHistos(selec_H.bare      ,                 rdfData       , gen = False )

  import pdb; pdb.set_trace()

  if newHb:
    selec_H_Bs           = createHistos(selec_H.bs          ,                 rdfBs         , gen = False,  variables = ["phiPi_m", "dsMu_m"])
    selec_H_B0           = createHistos(selec_H.b0          ,                 rdfB0         , gen = False,  variables = ["phiPi_m", "dsMu_m"])
    selec_H_Bplus        = createHistos(selec_H.bplus       ,                 rdfBplus      , gen = False,  variables = ["phiPi_m", "dsMu_m"])

    selec_H_Bs_2D        = create2DHistos(selec_H.bs        ,                 rdfBs         , gen = False )
    selec_H_B0_2D        = create2DHistos(selec_H.b0        ,                 rdfB0         , gen = False )
    selec_H_Bplus_2D     = create2DHistos(selec_H.bplus     ,                 rdfBplus      , gen = False )

  if constrained:
    #selec_C_Data_sf     = createHistos(baseline + score_cut + wrong_sign + low_mass,   rdfData        , gen = False, massPlot = True, variables = ["phiPi_m"])
    selec_H_Data_sf_kk   = createHistos(  baseline + score_cut + kk_wrong   ,      rdfData       , gen = False, variables = ["phiPi_m", "dsMu_m"])
    selec_H_Data_sf_pimu = createHistos(  baseline + score_cut + pimu_wrong ,      rdfData       , gen = False, variables = ["phiPi_m", "dsMu_m"])
    selec_H_Data_sf_2D   = create2DHistos(baseline + score_cut + wrong_sign ,      rdfData       , gen = False)
 
    #get the signflip scale by fitting the ds mass peak
  
    hbDummy = selec_H_Hb["phiPi_m"].Clone()
    if hb_scale_H == 0:
      hbDummy.Scale(0.0 )
    else:
      hbDummy.Scale(hb_scale_H / hbDummy.Integral())

    hRest    = selec_H_DsMu["phiPi_m"].Clone()
    hRest.Add( selec_H_DsStarMu["phiPi_m"].Clone())
    hRest.Add( selec_H_DsTau["phiPi_m"].Clone())
    hRest.Add( selec_H_DsStarTau["phiPi_m"].Clone())
    hRest.Add( hbDummy)

    global scale_H_kk;
    global scale_H_bkg;
    global scale_H_n;

    scale_H_kk, scale_H_bkg, scale_H_n = getSignflipRatio(selec_H_Data_sf_kk["phiPi_m"].Clone(),selec_H_Data_sf_pimu["phiPi_m"].Clone(),hRest,selec_H_Data["phiPi_m"].Clone())
 
  print("===> high mass region done...")



  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if args.pastNN: name += "_pastNN"
  else: name += ""



  for var in pastNN_2Dmodels.keys():

      myFile    = ROOT.TFile.Open(f"/work/pahwagne/RDsTools/fit/shapes/{var}_shapes_{name}_2D.root", "RECREATE")

      histos2D = {"dsTau":  selec_S_DsTau_2D[var],
               "dsStarTau": selec_S_DsStarTau_2D[var], 
               "dsMu":      selec_S_DsMu_2D[var], 
               "dsStarMu":  selec_S_DsStarMu_2D[var], 
               "hb"       : selec_S_Hb_2D[var], 
               "data"     : selec_S_Data_2D[var]}

      if newHb:
        histos2D["bs"]    = selec_S_Bs_2D[var] 
        histos2D["b0"]    = selec_S_B0_2D[var] 
        histos2D["bplus"] = selec_S_Bplus_2D[var] 

      if constrained:

        histos2D["data_sf"] = selec_S_Data_sf_2D[var] 
        toPass = histos2D.copy() 

        histos2DScaled = stacked2DPlot(toPass, var, hb_scale_S, mlow, mhigh, constrained = constrained, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        #normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, constrained = constrained, newHb = newHb)
        #normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, constrained = constrained, newHb = newHb, log = True)

      else:

        histos2D["combL"] = selec_S_DataL_2D[var] 
        histos2D["combR"] = selec_S_DataR_2D[var] 
        toPass = histos2D.copy() 

        histos2DScaled = stacked2DPlot(toPass, var, hb_scale_S, mlow, mhigh, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        #normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, newHb = newHb)
        #normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, log = True, newHb = newHb)

      myFile.WriteObject(histos2DScaled["dsMu"],      "dsMu"      )
      myFile.WriteObject(histos2DScaled["dsTau"],     "dsTau"     )
      myFile.WriteObject(histos2DScaled["dsStarMu"],  "dsStarMu"  )
      myFile.WriteObject(histos2DScaled["dsStarTau"], "dsStarTau" )
      myFile.WriteObject(histos2DScaled["hb"],        "hb"        )
      myFile.WriteObject(histos2DScaled["data"],      "data_obs"  ) #combine convention
      myFile.WriteObject(histos2DScaled["comb"],      "comb"      )

      writeDatacard(histos2DScaled, var, is_2D = True)

  for var in models.keys():

    #create root file for every variable which holds shapes
    myFile    = ROOT.TFile.Open(f"/work/pahwagne/RDsTools/fit/shapes/{var}_shapes_{name}.root", "RECREATE")
    myFile_1D = ROOT.TFile.Open(f"/work/pahwagne/RDsTools/fit/shapes/{var}_shapes_1D_{name}.root", "RECREATE")

    print(f"===> Producing stacked plot for variable: {var}") 
    if "gen" in var:
      #skip gen variables
      print("This is a gen variable... skip!")
      continue
 

    ######################################
    ## Plot variables except phiPi mass ##
    ######################################

    if (var != "phiPi_m") and (var!= "dsMu_m"):
  
      histos = {"dsTau":     selec_S_DsTau[var],
               "dsStarTau": selec_S_DsStarTau[var], 
               "dsMu":     selec_S_DsMu[var], 
               "dsStarMu":  selec_S_DsStarMu[var], 
               "hb"       : selec_S_Hb[var], 
               "data"     : selec_S_Data[var]}

      if newHb:
        histos["bs"]    = selec_S_Bs[var] 
        histos["b0"]    = selec_S_B0[var] 
        histos["bplus"] = selec_S_Bplus[var] 

      if constrained:

        #histos["data_sf"] = selec_S_Data_sf[var] 
        histos["data_sf_kk"]   = selec_S_Data_sf_kk[var] 
        histos["data_sf_pimu"] = selec_S_Data_sf_pimu[var] 
        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_S, scale_S_kk, scale_S_bkg, scale_S_n, mlow, mhigh, constrained = constrained, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_S_kk, constrained = constrained, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_S_kk, constrained = constrained, newHb = newHb, log = True)

      else:

        histos["combL"] = selec_S_DataL[var] 
        histos["combR"] = selec_S_DataR[var] 
        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_S, scale_S_kk, scale_S_bkg, scale_S_n, mlow, mhigh, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_S_kk, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_S_kk,constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, log = True, newHb = newHb)

      myFile.WriteObject(histosScaled["dsMu"],      "dsMu"      )
      myFile.WriteObject(histosScaled["dsTau"],     "dsTau"     )
      myFile.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile.WriteObject(histosScaled["hb"],        "hb"        )
      myFile.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      myFile.WriteObject(histosScaled["comb"],      "comb"      )

      myFile_1D.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile_1D.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile_1D.WriteObject(histosScaled["hb"],        "hb"        )
      myFile_1D.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      myFile_1D.WriteObject(histosScaled["comb"],      "comb"      )

      writeDatacard(histosScaled, var)
 
    if (var == "phiPi_m"): 
      histos = {"dsTau":    selec_C_DsTau[var],
               "dsStarTau": selec_C_DsStarTau[var], 
               "dsMu":      selec_C_DsMu[var], 
               "dsStarMu":  selec_C_DsStarMu[var], 
               "hb"       : selec_C_Hb[var], 
               "data"     : selec_C_Data[var]}

      if newHb:
        histos["bs"]    = selec_C_Bs[var] 
        histos["b0"]    = selec_C_B0[var] 
        histos["bplus"] = selec_C_Bplus[var] 

      if constrained:
 
        #histos["data_sf"]      = selec_C_Data_sf[var]
        histos["data_sf_kk"]   = selec_C_Data_sf_kk[var] 
        histos["data_sf_pimu"] = selec_C_Data_sf_pimu[var] 

        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_C, scale_C_kk, scale_C_bkg, scale_C_n, mlow, mhigh, constrained = constrained, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_C_kk,constrained = constrained)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_C_kk, constrained = constrained, log = True)

      else:

        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_C, scale_C_kk, scale_C_bkg, scale_C_n,mlow, mhigh, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_C_kk, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_C_kk, constrained = constrained, fakemass = fakemass, A = A, B = B, C = C, S = S, log = True)

      #stackedPlot(histos, var, hb_scale, fakemass,A,B,C,S, log = True)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S)
      #normPlot(dict(list(histos.items())[:-1]), var, fakemass, A,B,C,S, log = True)
      

      myFile.WriteObject(histosScaled["dsMu"],      "dsMu"      )
      myFile.WriteObject(histosScaled["dsTau"],     "dsTau"     )
      myFile.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile.WriteObject(histosScaled["hb"],        "hb"        )
      myFile.WriteObject(histosScaled["data"],      "data_obs"  ) #combine covention
      myFile.WriteObject(histosScaled["comb"],      "comb"      )

      myFile_1D.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile_1D.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile_1D.WriteObject(histosScaled["hb"],        "hb"        )
      myFile_1D.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      myFile_1D.WriteObject(histosScaled["comb"],      "comb"      )

      writeDatacard(histosScaled, var)

    if (var == "dsMu_m"): 
      histos = {"dsTau":    selec_H_DsTau[var],
               "dsStarTau": selec_H_DsStarTau[var], 
               "dsMu":      selec_H_DsMu[var], 
               "dsStarMu":  selec_H_DsStarMu[var], 
               "hb"       : selec_H_Hb[var], 
               "data"     : selec_H_Data[var]}

      if newHb:
        histos["bs"]    = selec_H_Bs[var] 
        histos["b0"]    = selec_H_B0[var] 
        histos["bplus"] = selec_H_Bplus[var] 

      if constrained:
 
        #histos["data_sf"]      = selec_H_Data_sf[var]
        histos["data_sf_kk"]   = selec_H_Data_sf_kk[var] 
        histos["data_sf_pimu"] = selec_H_Data_sf_pimu[var] 

        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_H, scale_H_kk, scale_H_bkg, scale_H_n, mlow_H, mhigh_H, constrained = constrained, bs = bs_in_hb_H, b0 = b0_in_hb_H, bplus = bplus_in_hb_H, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_H_kk, constrained = constrained)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_H_kk, constrained = constrained, log = True)

      else:

        toPass = histos.copy() 

        histosScaled = stackedPlot(toPass, var, hb_scale_H, scale_H_kk, scale_H_bkg, scale_H_n,mlow_H, mhigh_H, constrained = constrained, fakemass = fakemass, A = A_H, B = B_H, C = C_H, S = S_H, bs = bs_in_hb_H, b0 = b0_in_hb_H, bplus = bplus_in_hb_H, newHb = newHb)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_H_kk, constrained = constrained, fakemass = fakemass, A = A_H, B = B_H, C = C_H, S = S_H)
        normPlot({key: toPass[key] for key in histos.keys() if ((key != "data") and (key != "comb"))}, var, scale_H_kk, constrained = constrained, fakemass = fakemass, A = A_H, B = B_H, C = C_H, S = S_H, log = True)

      myFile.WriteObject(histosScaled["dsMu"],      "dsMu"      )
      myFile.WriteObject(histosScaled["dsTau"],     "dsTau"     )
      myFile.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile.WriteObject(histosScaled["hb"],        "hb"        )
      myFile.WriteObject(histosScaled["data"],      "data_obs"  ) #combine covention
      myFile.WriteObject(histosScaled["comb"],      "comb"      )

      myFile_1D.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      myFile_1D.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      myFile_1D.WriteObject(histosScaled["hb"],        "hb"        )
      myFile_1D.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      myFile_1D.WriteObject(histosScaled["comb"],      "comb"      )

      writeDatacard(histosScaled, var)

#createPlots(ma_cut, constrained = args.constrained, newHb = False)
createPlots(base_wout_tv, constrained = args.constrained, newHb = args.newHb)
#createPlots(base, constrained = args.constrained, newHb = args.newHb)

