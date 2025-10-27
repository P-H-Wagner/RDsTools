import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
ROOT.ROOT.EnableImplicitMT() 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio, fitAnotherVar
from helper import * 
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels, special_models, special_models_q2_coll
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
parser.add_argument("--fit",         required = True,     help = "Specify fitter 'cons' (constrained) or 'uncons' (unconstrained)") 
parser.add_argument("--hb",          required = True,     help = "Specify hb bkg modeling 'old' (inclusive) or 'new' (exclusive)") 
parser.add_argument("--nn",          required = True,     help = "Specify 'before' or 'past' neural network plots ") 
parser.add_argument("--hammer",      required = True,     help = "Specify 'true' or 'false' to apply hammer weights") 
parser.add_argument("--hammer_sys",  required = True,     help = "Specify 'true' or 'false' to save weight variation shapes") 
parser.add_argument("--prod",        required = True,     help = "Specify '24' or '25' to specify the data production") 
parser.add_argument("--trigger",     required = True,     help = "Specify mu7 or mu9 for the trigger") 
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
parser.add_argument("--control",                          help = "If given, run control plots, either 'highmass', 'sb'") 
parser.add_argument("--cut",                              help = "Cut on the discriminator score 5. ") 
parser.add_argument("--findcut",     action='store_true', help = "If given, we run thecut scan") 
args = parser.parse_args()

# python plotter.py -f cons -c old -n past -s true 


if args.fit not in ["cons", "uncons"]:
  raise ValueError ("Error: Not a valid key for --fit, please use 'cons' or 'uncons' (all lowercase!)")
else: constrained = (args.fit == "cons")

if args.hb not in  ["old", "new"]:
  raise ValueError ("Error: Not a valid key for --hb, please use 'old' or 'new' (all lowercase!)")
else: newHb       = (args.hb == "new")

if args.nn not in  ["before", "past"]:
  raise ValueError ("Error: Not a valid key for --nn, please use 'before' or 'past' (all lowercase!)")
else: pastNN      = (args.nn == "past")

if args.hammer not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --hammer, please use 'true' or 'false' (all lowercase!)")
else: hammer_central = (args.hammer == "true")

if args.hammer_sys not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --sys, please use 'true' or 'false' (all lowercase!)")
else: hammer_sys = (args.hammer_sys == "true")

if args.prod not in ["24", "25"]:
  raise ValueError ("Error: Not a valid key for --prod, please use '24' or '25'")
else: prod = args.prod

if args.trigger not in ["mu7", "mu9"]:
  raise ValueError ("Error: Not a valid key for --trigger, please use 'mu7' or 'mu9'")
else: trigger = args.trigger

if args.bdt not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt , please use 'true' or 'false' (all lowercase!)")
else: sf_weights = (args.bdt == "true")
print("sf_weights are:", sf_weights)

if args.cut  : score_cut = f" && (score5 <= {args.cut}) "
else         : score_cut = ""#"&& (dsMu_deltaR < 0.3) && (abs(cosPiK1) > 0.6) && (cosMuW_lhcb_alt > 0.6)"

print(score_cut)

if args.debug: debug = 50000
else         : debug = None



print(f"====> Running {args.fit} fit with {args.hb} Hb background, {args.nn} neural network on trigger {trigger}")
if pastNN: models.update(pastNN_models)

if hammer_sys:
  print(f"====> Adding the following systematics: {systematics_scalar} and {systematics_vector}") #defined in helper.py
  sys_scalar = systematics_scalar
  sys_vector = systematics_vector

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/20_10_2025_09_57_12/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)

if prod == "25":

  baseline      = base_wout_tv_25

  bdt_data   = bdt_data_25
  sig_cons   = sig_cons_25
  hb_cons    = hb_cons_25
  data_cons  = data_cons_25

  if trigger == "mu7": 
    cons_pastNN_25 = cons_pastNN_25_mu7 
  else:                
    cons_pastNN_25 = cons_pastNN_25_mu9

  sig_cons_pastNN     =  cons_pastNN_25["sig"]    
  hb_cons_pastNN      =  cons_pastNN_25["hb"]     
  data_cons_pastNN    =  cons_pastNN_25["data"]   

else:

  baseline      = base_wout_tv_24

  bdt_data   = bdt_data_24
  sig_cons   = sig_cons_24
  hb_cons    = hb_cons_24
  data_cons  = data_cons_24

  sig_cons_pastNN     =  sig_cons_pastNN_24    
  hb_cons_pastNN      =  hb_cons_pastNN_24     
  data_cons_pastNN    =  data_cons_pastNN_24   

skimmed = True
if args.control:
  if args.control == "highmass":
    skimmed = True
  elif args.control == "sb":
    skimmed = True
  else: 
    raise ValueError ("Error: Not a valid key for --control, please use 'highmass', 'sb' ")

control = args.control
print("Using baseline selection:", baseline)

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

# get date and time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S") 

#where to save plots and datacards
if args.findcut:
  toSave_plots        = f"/work/pahwagne/RDsTools/plots/cut_plots/{args.cut}/"
  toSave_cards        = f"/work/pahwagne/RDsTools/fit/cut_cards/{args.cut}/"
  shapes_folder       = f"/work/pahwagne/RDsTools/fit/cut_shapes/{args.cut}/"
  shapes_folder_blind = f"/work/pahwagne/RDsTools/fit/cut_shapes/{args.cut}/blind/"



else:
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

with open( toSave_plots + f"/info.txt", "a") as f:
  f.write( f"Using sf weights? {sf_weights} \n")

if trigger == "mu7":
  data_selec   = " && (mu7_ip4 == 1)"
  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
else:
  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "



# set the to be splitter variable and binning
#split   = "cosMuW_lhcb_alt"
split   = "q2_coll"
#split   = "q2_lhcb_alt"
#split   = "class"
#split   = "phiPi_m"
#split   = "e_star_reco_weighted"
#split = "score1"
#split = "score2"
#binning = [[0,0.03],[0.03,0.06],[0.06,0.1],[0.1,0.2],[0.2,0.4],[0.4,0.6],[0.6,1.0]] #score 3
#binning = [[0,0.02],[0.02,0.1],[0.1,0.4],[0.4,1]] #score 2
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,8.5],[8.5,9],[9,10],[10,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[9,12]]
#binning = [[-1,-0.5],[-0.5, -0.0], [0,0.5],[0.5,1]] #cosmuw
#binning = [[1.94, 1.95],[1.95, 1.96], [1.96, 1.98], [1.98,2.0]] #mass]
#binning = [[0,0.5], [0.5,1.0],[1.0, 1.5], [1.5, 2.0],[2.0, 3.0]] # estar
#binning = [[0,7],[7,12]]
binning = [[-100000,1000000]]
#binning = [[-0.01,0.01],[0.99,1.01],[1.99,2.01],[2.99,3.01],[3.99,4.01],[4.99,5.01]]
#binning = [[0.99,1.01]]
#binning = [[4.99,5.01]]
#binning = [[-99,99]]

################
## 2D binning ##
################

#in this case the splitter is irrelevant

#binning = [
#"(class == 0) && ((q2_coll < 7) || (q2_coll > 9) )",
#"(class == 1) && ((q2_coll < 7) || (q2_coll > 9) )",
#"(class == 2) && ((q2_coll < 7) || (q2_coll > 9) )",
#"(class == 3) && ((q2_coll < 7) || (q2_coll > 9) )",
#"(class == 4) && ((q2_coll < 7) || (q2_coll > 9) )",
#"(class == 5) && ((q2_coll < 7) || (q2_coll > 9) )",
#
#"(class == 0) && ((q2_coll > 7) && (q2_coll < 9) )",
#"(class == 1) && ((q2_coll > 7) && (q2_coll < 9) )",
#"(class == 2) && ((q2_coll > 7) && (q2_coll < 9) )",
#"(class == 3) && ((q2_coll > 7) && (q2_coll < 9) )",
#"(class == 4) && ((q2_coll > 7) && (q2_coll < 9) )",
#"(class == 5) && ((q2_coll > 7) && (q2_coll < 9) )",
#
#]

##############################
# Load chain into RDataFrame #
##############################

def getRdf(dateTimes, debug = None, skimmed = None, rdfSys = False, sf_weights = None, bph_part = None):

  print(dateTimes)
  chain = ROOT.TChain("tree")

  #if ((not pastNN) and (not isinstance(dateTimes, list))):
  #  print("dateTimes must be a list of strings")

  if (pastNN and not rdfSys):
    print(f"picking past NN files ...")
    #dateTimes here is a string like: "data_26Sep..._cons"
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/{dateTimes}/*"
    print(files)
    n = chain.Add(files)
    print(f"Chaining {n} files for this channel")

    if debug: 
      reduced_chain = chain.CloneTree(debug); 
      rdf = ROOT.RDataFrame(reduced_chain,debug)
      return (reduced_chain,rdf)

    else:     
      rdf = ROOT.RDataFrame(chain)
      return (chain,rdf)
 
  if (not pastNN and rdfSys):
    print(f"picking pre NN files with hammer weights...")
    #dateTimes here is a string like: "data_26Sep..._cons"
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/{dateTimes}/*"
    print(files)
    n = chain.Add(files)
    print(f"Chaining {n} files for this channel")

    if debug: 
      reduced_chain = chain.CloneTree(debug); 
      rdf = ROOT.RDataFrame(reduced_chain,debug)
      return (reduced_chain,rdf)

    else:     
      rdf = ROOT.RDataFrame(chain)
      return (chain,rdf)


  if sf_weights:
    files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data}/*" #test
    print("picking sf weighted files: ", files)
    chain.Add(files)
    if debug:
      reduced_chain = chain.CloneTree(debug); 
      rdf = ROOT.RDataFrame(reduced_chain, debug)
      return(reduced_chain, rdf)
    else:
      rdf = ROOT.RDataFrame(chain)
      return (chain,rdf)

  for dateTime in dateTimes:

    if bph_part is not None:
      print("Restrict plotting to bph part: ", bph_part)
      # if bph_part given, restrict to it.
      if dateTime != dateTimes[bph_part-1]: 
        print("skip this bph part ...")
        continue

    if skimmed:
      print("picking skimmed flatNano")
      files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/*" #skimmed files
      #files =  f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/{dateTime}/skimmed_bkg_{dateTime}.root"  # data skimmed with kkpimu > Bs for closure
      print(f"Appending {files}")


    else:
      #access the flat ntuples
      files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test
      print(f"Appending {files}")

    #chain them all
    chain.Add(files)

  if debug: 
    reduced_chain = chain.CloneTree(debug); 
    rdf = ROOT.RDataFrame(reduced_chain,debug)
    return (reduced_chain,rdf)

  else:     
    rdf = ROOT.RDataFrame(chain)
    print("rdf created")
    return (chain,rdf)

########################################
# Assign selections depending on the   #
# MC signal ID                         #
########################################

class selections:

  def __init__(self, selec):



    self.hb =        selec + mc_selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1 ) "
    self.bs =        selec + mc_selec + " && (300 <= gen_sig) && (gen_sig < 400) "
    self.b0 =        selec + mc_selec + " && (200 <= gen_sig) && (gen_sig < 300) "
    self.bplus =     selec + mc_selec + " && (100 <= gen_sig) && (gen_sig < 200) "
    self.dsMu =      selec + mc_selec + " && (gen_sig == 0)" 
    self.dsTau =     selec + mc_selec + " && (gen_sig == 1)"
    self.dsStarMu =  selec + mc_selec + " && (gen_sig == 10)" 
    self.dsStarTau = selec + mc_selec + " && (gen_sig == 11)"
    self.bare =      selec + data_selec

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
 
  print("Signal Region is:", signalRegion)

  return mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, signalRegion, anti_signalRegion, leftSB, rightSB

#########################################
## INPUT SAMPLES                       ##
#########################################

# Create rdf from tree
print(f" ===> Start creating RDataFrames")


if pastNN:

  print("Producing post-NN plots")

  #save the NN name
  import re
  pattern = r"\d{2}[A-Za-z]{3}\d{4}_\d{2}h\d{2}m\d{2}s"
  try:
    match = re.search(pattern,sig_cons_pastNN)
    nnModel = match.group(0)
  except: 
    nnModel = sig_cons_pastNN.split("_")[1]


  if constrained:

    sig_sample = sig_cons_pastNN
    print("sample is:", sig_cons_pastNN)
    print("constrained data after NN not processed yet..")

    chainSigSB, rdfSigSB     = getRdf(sig_cons_pastNN                 , debug = debug)
    chainSig,   rdfSig       = getRdf(sig_cons_pastNN                 , debug = debug)
    chainHb,    rdfHb        = getRdf(hb_cons_pastNN                  , debug = debug)
    _,          rdfHbHammer  = getRdf(hb_cons_pastNN                  , debug = debug)
    chainData,  rdfData      = getRdf(data_cons_pastNN                , debug = debug)

    print("rdf data has events: ", rdfData.Count().GetValue() )
    print("rdf sig has events: ", rdfSig.Count().GetValue() )

  else:

    if hammer_central: sig_sample = sig_unc_hammer
    else:              sig_sample = sig_unc_pastNN

    chainSigSB, rdfSigSB     = getRdf(sig_sample                       , rdfSys = hammer_central)#, debug = 10)
    chainSig,   rdfSig       = getRdf(sig_sample                       , rdfSys = hammer_central)#, debug = 10)
    chainHb,    rdfHb        = getRdf(hb_unc_pastNN                    )#, debug = 10)
    chainData,  rdfData      = getRdf(data_unc_pastNN                  )#, debug = 10)


else:

  print("Producing pre-NN plots")
  if constrained:

    if hammer_central: sig_cons = sig_cons_hammer_25

    chainSigSB, rdfSigSB     = getRdf(sig_cons                        , debug = debug, skimmed = skimmed, rdfSys = hammer_central )# , skimmed = baseline_name)#, debug = 1)
    chainSig,   rdfSig       = getRdf(sig_cons                        , debug = debug, skimmed = skimmed, rdfSys = hammer_central )# , skimmed = baseline_name)#, debug = 1)
    chainHb,    rdfHb        = getRdf(hb_cons                         , debug = debug, skimmed = skimmed )# , skimmed = baseline_name)#, debug = 1)
    chainData, rdfData       = getRdf(data_cons    ,sf_weights = sf_weights , debug = debug, skimmed = skimmed )# , skimmed = baseline_name)#, debug = 1)
    #print("---------------> rdf has events: ", rdfData.Count().GetValue() )


  else:
  
    chainSigSB, rdfSigSB     = getRdf(sig_unc      ,debug = 200000) #    , skimmed = baseline_name)#, debug = 10)
    chainSig,   rdfSig       = getRdf(sig_unc      ,debug = 200000) #    , skimmed = baseline_name)#, debug = 10)
    chainHb,    rdfHb        = getRdf(hb_unc       ,debug = 200000) #    , skimmed = baseline_name)#, debug = 10)
    chainData,  rdfData      = getRdf(data_unc     ,debug = 200000) #    , skimmed = baseline_name)#, debug = 10)
  
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
              "data_sf_both":    ROOT.kGray+1,
              "data_sf_one":    ROOT.kGray+1,
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
  nReco1 = chain.GetEntries(selection + f"&& ({name}_reco_1 > {start}) && ({name}_reco_1 < {stop}) && (gen_{name} > {start}) && (gen_{name} < {stop}) && (abs({name}_reco_1 - gen_{name}) < abs({name}_reco_2 - gen_{name})) && ({name}_reco_1 == {name}_reco_1)")
    
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
 
def addSystematics(hist_dict, var, selec_DsTau, selec_DsMu, selec_DsStarTau, selec_DsStarMu ):

   # if there are systematics, include them in the scaling
   for s in sys_scalar:
     for direc in ["Up", "Down"]:
       hist_dict["dsTau_"     + s + direc] = selec_DsTau    [var + "_" + s + direc]
       hist_dict["dsMu_"      + s + direc] = selec_DsMu     [var + "_" + s + direc]

   for s in sys_vector:
     for direc in ["Up", "Down"]:
       hist_dict["dsStarTau_" + s + direc] = selec_DsStarTau[var + "_" + s + direc]
       hist_dict["dsStarMu_"  + s + direc] = selec_DsStarMu [var + "_" + s + direc]

   return hist_dict


#########################################
## CREATE DEFAULT HISTOS               ##
#########################################

def createHistos(selection,rdf, linewidth = 2, gen = True, data = False , variables = None, hammer_central = False, hammer_sys = False, sig = None, sf_weights = None, region = None, massfit = False):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"
  "If <variables> are given, only these variables from <models> are created" #f.e. for the phiPi mass or DsMu with different selection

  print("we are here in bin", region)


  histos = {}

  if pastNN:           models.update(pastNN_models)
  if massfit == False: models.update(modelsSR)

  for var,model in models.items(): 

      print("var is", var, "split is", split)

      ################################
      # Special binnings             #
      ################################
      if var == "score1" and (split == "class") and region != None:
        
        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")
  
        model = special_models[var + f"_bin{region}" ]

      if var == "score1" and (split == "q2_coll" or split == "q2_lhcb_alt") and region != None:
        
        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")
  
        model = special_models_q2_coll[var + f"_bin{region}" ]                                     


      ##############################################
      # Fill only variables, if explicitly given!  #
      ##############################################

      if variables:
        if var not in variables: 
          print(f"cannot find {var}");     
          continue 

      ##############################################
      # Skip gen variables ?                       #
      ##############################################

      if "gen" in var and not gen:
        #skip gen variables
        print("This is a gen variable... skip!")
        continue

      tofill = var

      ##############################################
      # Shift mass peak of MC by half permill      #
      ##############################################

      if (data == False and var == "phiPi_m"):
        tofill = f"(run==1) * 0.9995 * phiPi_m + (run!=1) * phiPi_m"
        print("correcting and plotting phiPi mass..")

        if hammer_central:
          #this only happens for signal
          #central weights

          print(averages)
          central_av = averages[ "central_w_" + sig]
           
          histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], "m_corr","central_w") 
          histos[var].Scale(1.0 / central_av) 
 

        elif sf_weights:
          # this only happens for data
          # sf_weight is a smart string
          histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], "m_corr" ,"sf_weights")

        else:
          #fill without any weights
          histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], "m_corr" )


        if hammer_sys:
          #add also systematical shape variations as variables 

          if "star" not in sig: sys_dir = sys_scalar #only take e1 - e6 for scalar signals (BCL)
          else:                 sys_dir = sys_vector #take e1-e10

          for s in sys_dir:

            s_up   = s + "_up"
            s_down = s + "_down"

            var_up_av   = averages[s_up   + "_" + sig]
            var_down_av = averages[s_down + "_" + sig]

            func_up   = s_up   + f" / ({ var_up_av   })" 
            func_down = s_down + f" / ({ var_down_av })" 

            # fill histogram with weight "s"
            histos[var + "_" + s + "Up"]   = rdf.Filter(selection).Define("m_corr", tofill).Histo1D(model[0], "m_corr", s + "_up"   )
            histos[var + "_" + s + "Up"].Scale(1.0 / var_up_av)

            histos[var + "_" + s + "Down"]   = rdf.Filter(selection).Define("m_corr", tofill).Histo1D(model[0], "m_corr", s + "_down"   )
            histos[var + "_" + s + "Down" ].Scale(1.0 / var_down_av)


      ##############################################
      # Histo creation of other variables          #
      ##############################################
      else:
        print("Filling for variable ", var)

        ##############################################
        # Histo creation of systematic up/down       #
        ##############################################

        if hammer_central: 

          central_av = averages[ "central_w_" + sig]
            
          histos[var] = rdf.Filter(selection).Histo1D(model[0], tofill,"central_w") 
          histos[var].Scale(1.0 / central_av) 
 
        elif (sf_weights):
          # this only happens for data
          # sf_weight is a smart string
          histos[var] = rdf.Filter(selection).Define("m_corr",tofill).Histo1D(model[0], tofill, "sf_weights")


        else:
          histos[var] = rdf.Filter(selection).Histo1D(model[0], tofill)

        if hammer_sys:
          #add also systematical shape variations as variables 
          # fill sysmetatic up and down variations
          if "star" not in sig: sys_dir = sys_scalar     #only take e1 - e6 for scalar signals (BCL)
          else:                 sys_dir = sys_vector     #take e1-e10
          for s in sys_dir:
 
            s_up   = s + "_up"
            s_down = s + "_down"

            var_up_av   = averages[s_up   + "_" + sig]
            var_down_av = averages[s_down + "_" + sig]

            func_up   = s_up   + f" / ({ var_up_av   })" 
            func_down = s_down + f" / ({ var_down_av })" 

            # fill histogram with weight "s"
            histos[var + "_" + s + "Up"]   = rdf.Filter(selection).Histo1D(model[0], tofill, s + "_up"   )
            histos[var + "_" + s + "Up"].Scale(1.0 / var_up_av)
            histos[var + "_" + s + "Down"]   = rdf.Filter(selection).Histo1D(model[0], tofill, s + "_down"   )
            histos[var + "_" + s + "Down" ].Scale(1.0 / var_down_av)


      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetBinContent(histos[var].GetMaximumBin()))
      
      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(linewidth)

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



def stackedPlot(histos2, var, hb_scale, abc = None, scale_kk = None, scale_pimu = None, scale_bkg = None,  scale_rest = None, scale_n = None, log = False, fakemass = None, A = None, B = None, C = None, S = None, bs = None, b0 = None, bplus = None, region = None, blind = False):

  histos = {key: hist.Clone() for key, hist in histos2.items()}
  

  #print("==========> number of tau signals right after calling stackedPlot", histos["dsTau"].Integral())
  #print("==========> number of tau* signals right after calling stackedPlot", histos["dsStarTau"].Integral())
  #print("==========> number of mu signals right after calling stackedPlot", histos["dsMu"].Integral())
  #print("==========> number of mu* signals right after calling stackedPlot", histos["dsStarMu"].Integral())
  #print("==========> number of pimu  wrong right after calling stackedPlot", histos["data_sf_pimu"].Integral())
  #print("==========> number of kk wrong right after calling stackedPlot", histos["data_sf_kk"].Integral())
  #print("==========> number of data events right calling stackedPlot", histos["data"].Integral())

  color, labels = getColorAndLabelSignalDistinction("stacked")

  #supress data plotting in SR. region is a string of type: q2_coll_ch0
  dataBlind = False
  if region[-1] in ["5","6","7"]:
    dataBlind = True 



  for i,key in enumerate(histos.keys()):

    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillColor(color[key])
    histos[key].SetLineColor(color[key])
    if key == "data":
       histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["dsMu"].GetNbinsX()

  # stacked histo
  hs = ROOT.THStack(var,"")
  # error histo
  hErr = ROOT.TH1D(histos["dsMu"])  
  hErr.Reset("ICES") 

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
  hErr.SetFillStyle(3144)
  hErr.SetMarkerStyle(0)

  #scale hb to other mc

  # if there is no hb in this bin, there is nothing to scale!
  if histos["hb"].Integral() != 0: 
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
  print("==========> number of tau* signals before scaling", histos["dsStarTau"].Integral())
  print("==========> number of mu* signals before scaling", histos["dsStarMu"].Integral())


  ###################################
  ## Combinatorial treatment       ##
  ###################################

  if constrained:

    # normalize
    hComb      = histos["data_sf_kk"].Clone()
    hComb.Scale(scale_kk)

    hComb_kk   = hComb.Clone()
    hComb_kk.SetLineColor(ROOT.kCyan +4)
    hComb_kk.SetFillColor(ROOT.kCyan +4)
    hComb_kk.SetFillStyle(3444)

 
    hComb_pimu = histos["data_sf_pimu"].Clone()
    hComb_pimu.Scale(scale_pimu)

    hComb_pimu.SetLineColor(ROOT.kMagenta -3)
    hComb_pimu.SetFillStyle(3244)
    hComb_pimu.SetFillColor(ROOT.kMagenta -3)

    hComb.Add(hComb_pimu)

    #pdb.set_trace()

    #hComb = histos["data_sf_one"].Clone()
    #hComb.Scale(scale_bkg)

    print("hcomb has now integral: ", hComb.Integral())


  else:     
    A = abc[0]
    B = abc[1]
    C = abc[2]
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
      for i in range(21):
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
        if nSig > 0:
          # if we have any signal or hb, normalize the proportions
          # if nSig = 0, there is nothing to normalize, all are zero
          #histos[key].Scale((histos["data"].Integral() - nComb) /nSig )
          histos[key].Scale(scale_rest) # NEW

    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    #print("nsig has integral:", nSig)

    for key in histos.keys(): # NEW
      #normalize to data!
      if ('data' not in key):
        histos[key].Scale( histos["data"].Integral() / (nSig + nComb))
        #print("scale_n inside :", histos["data"].Integral())

  else:
    #For sideband method, combinatorial is already normalized to data!
    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key):
        histos[key].Scale((histos["data"].Integral() - nComb) / nSig)

  #print("==========> number of tau* signals after scaling MC to data", histos["dsStarTau"].Integral())
  #print("==========> number of mu* signals after scaling MC to data", histos["dsStarMu"].Integral())
  #print("==========> number of pimu  wrong right after scaling mc to data ", histos["data_sf_pimu"].Integral())
  #print("==========> number of kk wrong right after scalinf mc to data", histos["data_sf_kk"].Integral())
  #print("==========> number of data avents after calling stackedPlot", histos["data"].Integral())

  # add comb to stack and err 
  hErr.Add(hComb)
  hs.Add(hComb)


  if newHb:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key) and ("Up" not in key) and ("Down" not in key):
        hErr.Add(histos[key])
        hs.Add(histos[key])
  else:
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key) and ("Up" not in key) and ("Down" not in key):
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
  leg.AddEntry(histos["comb"]                           ,'Comb. + Fakes'  ,'F' )
  leg.AddEntry(hErr                            ,'Stat. Uncer.'  ,'F' )
  
  #plot mainpad
  if log:
    #ROOT.gPad.SetLogy()
    #histos["data"].GetYaxis().SetTitle('Events')
    #histos["data"].GetYaxis().SetRangeUser(1e-3, hs.GetMaximum()*1000)
    #histos["data"].SetMinimum(1.001) # avoid idsplaying tons of comb bkg

    ROOT.gPad.SetLogx()

    yAxisTitle = "events"
    histos["data"].GetYaxis().SetTitle(yAxisTitle)
    histos["data"].GetYaxis().SetRangeUser(1e-3, histos["data"].GetBinContent(histos["data"].GetMaximumBin())*1.5)

  else:

    yAxisTitle = "events"
    histos["data"].GetYaxis().SetTitle(yAxisTitle)
    histos["data"].GetYaxis().SetRangeUser(1e-3, histos["data"].GetBinContent(histos["data"].GetMaximumBin())*1.5)

  #print("finally here , data : ", histos["data"].Integral() )
 
  #print("==========> number of tau signals inside stacked", histos["dsStarTau"].Integral())
  #for i in range(0, histos["data"].GetNbinsX() + 2):  # ROOT bins start from 1
  #        bin_center = histos["data"].GetBinCenter(i)
  #        bin_content = histos["data"].GetBinContent(i)
  #        print(f"Bin {i} (x={bin_center}): {bin_content}")

  #if we are in the signal region and do not blind the MC, then we blind the data!

  if dataBlind and not blind:
    #only draws axis, but not the data! :D
    histos["data"].Draw("AXIS")

  else:
    #draw data with uncertainty
    histos["data"].Draw("EP")

  #draw stacked histo
  hs.Draw("HIST SAME")

  #if constrained:
  #  #draw the two comb contributoins #new!!
  #  hComb_kk.Draw("HIST SAME")
  #  hComb_pimu.Draw("HIST SAME")

  #draw data again (used for plotting uses)
  if dataBlind and not blind:
    histos["data"].Draw("AXIS SAME")
  else:
    histos["data"].Draw("EP SAME")


  #draw stat uncertainty of MC, (can not be directly done with a stacked histo)
  hErr.Draw("E2 SAME")
  #hErr.Draw("E2 SAME")

  leg.Draw("SAME")
  ROOT.gPad.RedrawAxis()
  
  if trigger == "mu7": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  if trigger == "mu9": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = '20.7 fb^{-1}')
 
 
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
      if ('comb' not in key) and ('data' not in key) and ("hb" not in key) and ("Up" not in key) and ("Down" not in key):
        h = histos[key].Clone()
        h.Divide(hErr)
        norm_stack.Add(h)

  else: 
    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key) and ('b0' not in key) and ('bs' not in key) and ('bplus' not in key) and ("Up" not in key) and ("Down" not in key):
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

  if log:
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



  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  if region:
    name += f"_{region}" 

  if blind: 
    toSave = toSave_plots + "blind/"
  else: 
    toSave = toSave_plots

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  # save nn info
  with open( toSave + f"/info.txt", "a") as f:
    f.write(f"Stacked plot for variable: {var}")
    if pastNN: 
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    f.write( f" Baseline selection is: \n {baseline} \n")
    if region: f.write( f" Splitting variable {split} into regions {binning} \n")

    for key in histos.keys():
      f.write( f" Overall number of {key} entries: {histos[key].Integral()} \n")

      n_bins = histos[key].GetNbinsX() 
      for i in range(1, n_bins + 1):  # ROOT bins start from 1
        f.write( f" number of {key} entries in bin {i}: {histos[key].GetBinContent(i)} \n")


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


  # last but not least calculate the weights between sf and data
  # this is only relevant for the first ds+mu pt stacked histo
  weights = []
  n_bins = hComb.GetNbinsX() 
  for i in range(1, n_bins + 1):  # ROOT bins start from 1

    bin_low_edge  = hComb.GetBinLowEdge(i)
    bin_high_edge = hComb.GetBinLowEdge(i + 1)

    bin_comb = hComb.GetBinContent(i)
    bin_data = histos["data"].GetBinContent(i)
    if bin_comb > 0:
      weight   = bin_data / bin_comb
    else:
      weight = 1
    sub_str  = f" (({var} >= {bin_low_edge}) && ({var} < {bin_high_edge} )) * {weight}" 
    weights.append(sub_str)

  weight_str = " + ".join(weights)
  #print(weight_str)

  return returnHisto, weight_str

def shapesPlot(histos2, var, hb_scale, abc = None, scale_kk = None, scale_pimu = None, scale_bkg = None, scale_rest = None, scale_n = None, log = False, fakemass = None, A = None, B = None, C = None, S = None, bs = None, b0 = None, bplus = None, region = None, blind = False):

  histos = {key: hist.Clone() for key, hist in histos2.items()}
  color, labels = getColorAndLabelSignalDistinction("stacked")

  for i,key in enumerate(histos.keys()):
    print("my key in histos is: ", key)
    if ("Up" in key) or ("Down" in key): continue

    try: histos[key] = histos[key].GetValue()
    except: histos[key] = histos[key]

    histos[key].SetFillStyle(0)
    histos[key].SetLineColor(color[key])
    histos[key].SetLineWidth(3)

    if key == "data":
       histos[key].SetMarkerStyle(8) 

  # take as model
  bins  = histos["dsMu"].GetNbinsX()
  start = histos["dsMu"].GetXaxis().GetBinLowEdge(1)   # first bin
  stop  = histos["dsMu"].GetXaxis().GetBinUpEdge(bins) # last bin

  # mass histo for fakemass

  hComb = ROOT.TH1D(var,"",bins,start,stop)
  hComb.SetName("mass")
  hComb.GetXaxis().SetTitleOffset(1.3)
  hComb.GetYaxis().SetTitleSize(1*hComb.GetYaxis().GetTitleSize())
  hComb.GetXaxis().SetTitleSize(1*hComb.GetXaxis().GetTitleSize())  
  hComb.SetFillColor(ROOT.kGray)
  hComb.SetLineColor(ROOT.kGray)

  #scale hb to other mc
  print(hb_scale)
  if histos["hb"].Integral() != 0: 
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

    hComb = histos["data_sf_one"].Clone()
    hComb.Scale(scale_bkg)

    print("hcomb has now integral: ", hComb.Integral())

  else: 
    #do sideband method
    A = abc[0]
    B = abc[1]
    C = abc[2]

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
        if nSig > 0:
          # if we have any signal or hb, normalize the proportions
          # if nSig = 0, there is nothing to normalize, all are zero
          #histos[key].Scale((histos["data"].Integral() - nComb) /nSig )
          histos[key].Scale(scale_rest) # NEW

    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    #print("nsig has integral:", nSig)

    for key in histos.keys(): # NEW
      #normalize to data!
      if ('data' not in key):
        histos[key].Scale( histos["data"].Integral() / (nSig + nComb))
        #print("scale_n inside :", histos["data"].Integral())

  else:
    #For sideband method, combinatorial is already normalized to data!
    nSig = histos["dsMu"].Integral() + histos["dsTau"].Integral() + histos["dsStarMu"].Integral() + histos["dsStarTau"].Integral() + histos["hb"].Integral()

    for key in histos.keys():
      if ('comb' not in key) and ('data' not in key):
        histos[key].Scale((histos["data"].Integral() - nComb) / nSig)

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
  leg.AddEntry(histos["dsStarMu"]   ,'B_{s}#rightarrow D*_{s}#mu#nu' ,'L' )
  leg.AddEntry(histos["dsMu"]       ,'B_{s}#rightarrow D_{s}#mu#nu'  ,'L' )
  leg.AddEntry(histos["dsStarTau"]  ,'B_{s}#rightarrow D*_{s}#tau#nu','L' )
  leg.AddEntry(histos["dsTau"]      ,'B_{s}#rightarrow D_{s}#tau#nu' ,'L' )

  if newHb:
    leg.AddEntry(histos["bs"]       ,'B_{s}#rightarrow D_{s} + #mu ' ,'L' )
    leg.AddEntry(histos["b0"]       ,'B_{0}#rightarrow D_{s} + #mu ' ,'L' )
    leg.AddEntry(histos["bplus"]    ,'B_{+}#rightarrow D_{s} + #mu ' ,'L' )
  else:
    leg.AddEntry(histos["hb"]         ,'H_{b}#rightarrow D_{s} + #mu ' ,'L' )
  #leg.AddEntry(histos["data"]       ,'Data','LEP')
  leg.AddEntry(histos["comb"]                           ,'Comb. + Fakes'  ,'L' )
  
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

    if "data" in key: continue
    print(f"drawing variable {key}")

    #if ("comb" in key) or ("hb" in key and newHb == True): 
    #  #comes later, dont draw sidebands
    #  continue

    if i == 0:
      histos[key].Draw("HIST")
    else:
      histos[key].Draw("HIST SAME")

  leg.Draw("SAME")

  #CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  
  #saving
  c1.Modified()
  c1.Update()

  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  if region:
    name += f"_{region}" 

  toSave = toSave_plots + "/shapes/"

  if blind: 
    toSave += "blind/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  with open( toSave + f"/info.txt", "a") as f:
    # save nn info
    if pastNN: 
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    # write selection info
    f.write( f" Baseline selection is: \n {baseline} \n")
    if region: f.write( f" Splitting variable {split} into regions {binning} \n")
    # save the bin entry numbers in txt file
    for key in histos.keys():
      f.write( f" Overall number of {key} entries: {histos[key].Integral()} \n")

      n_bins = histos[key].GetNbinsX() 
      for i in range(1, n_bins + 1):  # ROOT bins start from 1
        f.write( f" number of {key} entries in bin {i}: {histos[key].GetBinContent(i)} \n")
      

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

    if key == "data":
       histos[key].SetMarkerStyle(8) 


  # mass histo for fakemass
  histos["data_sf_kk"]  .SetLineColor(ROOT.kGreen)
  histos["data_sf_pimu"].SetLineColor(ROOT.kCyan)
  histos["data_sf_both"].SetLineColor(ROOT.kBlue)
  histos["data_sf_one"].SetLineColor(ROOT.kOrange)

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

  leg.AddEntry(histos["data"]           ,'Right sign data'  ,'L' )
  leg.AddEntry(histos["data_sf_kk"]     ,'kk flip'          ,'L' )
  leg.AddEntry(histos["data_sf_pimu"]   ,'pimu flip'        ,'L' )
  leg.AddEntry(histos["data_sf_both"]   ,'both AND'         ,'L' )
  leg.AddEntry(histos["data_sf_one"]    ,'both OR'          ,'L' )
  
  #plot mainpad

  #normalize -> drawing with NORM doesnt work!!!
  histos["data"]        .Scale(1./histos["data"]        .Integral())  
  histos["data_sf_kk"]  .Scale(1./histos["data_sf_kk"]  .Integral())   
  histos["data_sf_pimu"].Scale(1./histos["data_sf_pimu"].Integral())  
  histos["data_sf_both"].Scale(1./histos["data_sf_both"].Integral())
  histos["data_sf_one"].Scale(1./histos["data_sf_one"].Integral())

  histos["data"].GetYaxis().SetTitle('a.u.')
  histos["data"].GetYaxis().SetRangeUser(1e-3, histos["data"].GetMaximum()*2.0)
 
  histos["data"]        .Draw("EP")
  histos["data_sf_kk"]  .Draw("HIST SAME") 
  histos["data_sf_pimu"].Draw("HIST SAME")
  histos["data_sf_both"].Draw("HIST SAME")
  histos["data_sf_one"].Draw("HIST SAME")

  leg.Draw("SAME")

  #saving
  c1.Modified()
  c1.Update()

  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  if region:
    name += f"_{region}" 

  toSave = toSave_plots + "/shapes_bkg/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  with open( toSave + f"/info.txt", "a") as f:
    # save nn info
    if pastNN: 
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    # write selection info
    f.write( f" Baseline selection is: \n {baseline} \n")
    if region: f.write( f" Splitting variable {split} into regions {binning} \n")
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
  
  if trigger == "mu7": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = '6.4 fb^{-1}')
  if trigger == "mu9": CMS_lumi(main_pad, 4, 0, cmsText = '     CMS', extraText = '        Preliminary', lumi_13TeV = '20.7 fb^{-1}')

  #saving
  c1.Modified()
  c1.Update()


  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if pastNN: 
    name += "_pastNN"
    folder = nnModel
  else: 
    name += ""
    folder = "/"

  toSave = f"/work/pahwagne/RDsTools/plots/normplots/{dt}/"

  if not os.path.exists(toSave): 
    os.makedirs(toSave)
    os.makedirs(toSave + "/log")  

  # save nn info
  with open( toSave + f"/info.txt", "a") as f:
    if pastNN:
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    f.write( f" Baseline selection is: \n {baseline} \n")


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

  # save nn info
  with open( toSave + f"/info.txt", "a") as f:

    if pastNN: 
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    f.write( f" Baseline selection is: \n {baseline} \n")

  canvas.SaveAs( toSave + f"/{name}.pdf")  
  canvas.SaveAs( toSave + f"/{name}.png")  
  print(f"DONE")


def writeShapes(hist_dict, outputFile, binned = False, channel = "placeholder :)"):

      #print(hist_dict["dsMu"])
      #print(outputFile)
      outputFile.WriteObject(hist_dict["dsMu"],      "dsMu"      + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["dsTau"],     "dsTau"     + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["dsStarMu"],  "dsStarMu"  + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["dsStarTau"], "dsStarTau" + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["hb"],        "hb"        + binned * f"_ch{channel}" )
      outputFile.WriteObject(hist_dict["data"],      "data_obs"  + binned * f"_ch{channel}" ) #combine convention
      outputFile.WriteObject(hist_dict["comb"],      "comb"      + binned * f"_ch{channel}" )

      if hammer_sys:
        # if there are systematics, include them in the scaling
        for s in sys_scalar:
          for direc in ["Up", "Down"]:
            outputFile.WriteObject(hist_dict["dsTau_"        + s + direc],"dsTau_"        + s + scalar_model + direc + binned * f"_ch{channel}"  )
            outputFile.WriteObject(hist_dict["dsMu_"         + s + direc],"dsMu_"         + s + scalar_model + direc + binned * f"_ch{channel}"  )

        for s in sys_vector:
          for direc in ["Up", "Down"]:
            outputFile.WriteObject(hist_dict["dsStarTau_"    + s + direc],"dsStarTau_"    + s + vector_model + direc + binned * f"_ch{channel}"  )
            outputFile.WriteObject(hist_dict["dsStarMu_"     + s + direc],"dsStarMu_"     + s + vector_model + direc + binned * f"_ch{channel}"  )

      #myFile_1D.WriteObject(histosScaled["dsStarMu"],  "dsStarMu"  )
      #myFile_1D.WriteObject(histosScaled["dsStarTau"], "dsStarTau" )
      #myFile_1D.WriteObject(histosScaled["hb"],        "hb"        )
      #myFile_1D.WriteObject(histosScaled["data"],      "data_obs"  ) #combine convention
      #myFile_1D.WriteObject(histosScaled["comb"],      "comb"      )


def writeBinnedDatacard(histos, var, region, splitter, digits = 5, blind = False):

  if constrained: name = "constrained"
  else: name = "unconstrained"

  if newHb: name += "_newHb"
  else: name += ""

  if pastNN: name += "_pastNN"
  else: name += ""

  if blind:
    toSave = toSave_cards + "/blind"
  else:
    toSave = toSave_cards

  if not os.path.exists(toSave): 
    os.makedirs(toSave)

  with open( toSave_plots + f"/info.txt", "a") as f:
    if pastNN:
      f.write( f" These plots use the following NN model: {nnModel} \n")
      f.write( f" And the following score cut: {score_cut} \n")
    f.write( f" Baseline selection is: \n {baseline} \n")
    f.write( f" Splitting variable {split} into regions {binning} \n")

  if hammer_sys:
    if args.findcut:
      temp = open("/work/pahwagne/RDsTools/fit/datacardTemplateSystematics_binned_cuts.txt", "rt")
    else:
      temp = open("/work/pahwagne/RDsTools/fit/datacardTemplateSystematics_binned.txt", "rt")

  else:
    if args.findcut:
      temp = open("/work/pahwagne/RDsTools/fit/datacardTemplate_binned_cuts.txt", "rt")
    else:
      temp = open("/work/pahwagne/RDsTools/fit/datacardTemplate_binned.txt", "rt")

  card_dir = toSave + f"/datacard_binned_{var}_{name}_{splitter}_ch{region}.txt"
  card = open(card_dir, "wt")

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
    if "HOOK_REGION"            in line: line = line.replace("HOOK_REGION",    str(region))
    if "HOOK_SPLITTER"          in line: line = line.replace("HOOK_SPLITTER",  splitter   )
    if "HOOK_RATES"             in line: line = line.replace("HOOK_RATES",     rates      )
    if "HOOK_DATA_RATE"         in line: line = line.replace("HOOK_DATA_RATE", dataStr    )
    if "HOOK_VAR"               in line: line = line.replace("HOOK_VAR",       var        )
    if "HOOK_NAME"              in line: line = line.replace("HOOK_NAME",      name       )
    if "HOOK_DATETIME"          in line: line = line.replace("HOOK_DATETIME",  dt )
    if "HOOK_CUT"               in line: line = line.replace("HOOK_CUT",       args.cut )
    if "HOOK_BLIND"             in line: line = line.replace("HOOK_BLIND",     blind * "blind" )
    card.write(line)

  temp.close()
  card.close()

  #return the path to the card
  return card_dir


def createBinnedPlots(splitter, regions, controlPlotsHighMass = None, controlPlotsSidebands = None):
  

  # splitter is the varible you want to use to define regions. F.e. q^2
  # regions is a list of intervals in the splitter variable, f.e. [[0,4],[4,8],[8,12]]
  # we will plot variables defined in models


  # this dictionary will hold a list of datacard paths (corresponding to the different regions) for every variable
  # since we define the regions in the variable {splitter} itself, we dont include the splitter variable to plot
  #cards       = {var: [] for var in models.keys() if ((var != splitter)  and (var != "phiPi_m"))}
  cards       = {var: [] for var in models.keys()}
  #cards_blind = {var: [] for var in models.keys() if ((var != splitter)  and (var != "phiPi_m"))}
  cards_blind = {var: [] for var in models.keys()}



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

  ################################################
  # Extract sigma from fit to MC DsMu signal to  #
  # define signal region.                        #
  ################################################

  # we freeze the signal region once and for all. 
  #selection_sigma = baseline + right_sign + low_mass + "&& (gen_sig == 0)" 

  # if we want to do control plots and thus our baseline is killing the signal, we just eyeball
  #if controlPlotsHighMass or controlPlotsSidebands: 
  #  sigma = 0.009 

  #else:
  #  sigma,   h            = getSigma(rdfSigSB, "phiPi_m", selection_sigma, dt)

  #also freeze sigma here!
  sigma = 0.009

  # define signal region and sidebands :)
  mlow,   mhigh,   mlow2,   mhigh2,   mlow3,   mhigh3,   signalRegion,   anti_signalRegion,   leftSB,   rightSB   = getRegions(sigma)

  ################################
  # Fit the Ds mass for signflip #
  ################################

  # We fit the Ds mass on the most general region (just apply the baseline selection) 
  # and get form it the ratio of yields.
  # Later, if we continue with the binning loop, we can just apply this ratio, the selection
  # is automatically applied since the to-be-plotted histograms already care the selection  in its belly


  if constrained:
 
    selection_massfit = baseline + right_sign 
    selec_massfit     = selections(selection_massfit)
    #hb_scale_massfit  = getHbScale(selec_massfit.hb , selec_massfit.dsMu )

    ## Signal and Hb 
    #selec_M_DsMu            = createHistos(selec_massfit.dsMu        ,    rdfSig       , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsmu"     )
    #selec_M_DsMu_woHammer   = createHistos(selec_massfit.dsMu        ,    rdfSig       , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False)
    #selec_M_DsTau           = createHistos(selec_massfit.dsTau       ,    rdfSig       , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dstau"    )
    #selec_M_DsStarMu        = createHistos(selec_massfit.dsStarMu    ,    rdfSig       , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsstarmu" )
    #selec_M_DsStarTau       = createHistos(selec_massfit.dsStarTau   ,    rdfSig       , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsstartau")
    #selec_M_Hb              = createHistos(selec_massfit.hb          ,    rdfHb        , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False)
    #selec_M_Mu_in_Hb        = createHistos(selec_massfit.dsMu        ,    rdfHb        , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False)
    #selec_M_Data            = createHistos(selec_massfit.bare        ,    rdfData      , variables = ["phiPi_m", "dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False)

    selec_M_DsMu            = createHistos(selec_massfit.dsMu        ,    rdfSig       , variables = ["phiPi_m"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsmu"     , massfit = True)
    selec_M_DsMu_woHammer   = createHistos(selec_massfit.dsMu        ,    rdfSig       , variables = ["phiPi_m"] , gen = False                                                                            , massfit = True)
    selec_M_DsTau           = createHistos(selec_massfit.dsTau       ,    rdfSig       , variables = ["phiPi_m"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dstau"    , massfit = True)
    selec_M_DsStarMu        = createHistos(selec_massfit.dsStarMu    ,    rdfSig       , variables = ["phiPi_m"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsstarmu" , massfit = True)
    selec_M_DsStarTau       = createHistos(selec_massfit.dsStarTau   ,    rdfSig       , variables = ["phiPi_m"] , gen = False, hammer_central = hammer_central, hammer_sys = hammer_sys,sig = "dsstartau", massfit = True)
    selec_M_Hb              = createHistos(selec_massfit.hb          ,    rdfHb        , variables = ["phiPi_m"] , gen = False                                                                            , massfit = True)
    selec_M_Mu_in_Hb        = createHistos(selec_massfit.dsMu        ,    rdfHb        , variables = ["phiPi_m"] , gen = False                                                                            , massfit = True)
    selec_M_Data            = createHistos(selec_massfit.bare        ,    rdfData      , variables = ["phiPi_m"] , gen = False, data = True                                                               , massfit = True)

    selec_M_DsTau_blind     = { key: selec_M_DsTau[key].Clone()      for key in selec_M_DsTau.keys()     }
    selec_M_DsStarTau_blind = { key: selec_M_DsStarTau[key].Clone()  for key in selec_M_DsStarTau.keys() }

    #selec_M_DsTau["phiPi_m"].Scale(0.3)
    #selec_M_DsStarTau["phiPi_m"].Scale(0.3)


    #since the dsmu signal in Hb is not hammered, we calc. the scale for unhammered dsmu and hb. The effects should cancel! 
    hb_ratio_massfit = selec_M_Hb["phiPi_m"].Integral() / selec_M_Mu_in_Hb["phiPi_m"].Integral()
    hb_scale_massfit = selec_M_DsMu_woHammer["phiPi_m"].Integral() * hb_ratio_massfit 
     
 
    for key in selec_M_DsTau_blind.keys()    : selec_M_DsTau_blind[key]    .Scale(blind_scalar)    
    for key in selec_M_DsStarTau_blind.keys(): selec_M_DsStarTau_blind[key].Scale(blind_vector)    

    selec_M_Hb_new = {}  

    if newHb:
      selec_M_Bs           = createHistos(selec_massfit.bs           ,    rdfBs         , variables = ["phiPi_m"] , gen = False, massfit = True)
      selec_M_B0           = createHistos(selec_massfit.b0           ,    rdfB0         , variables = ["phiPi_m"] , gen = False, massfit = True)
      selec_M_Bplus        = createHistos(selec_massfit.bplus        ,    rdfBplus      , variables = ["phiPi_m"] , gen = False, massfit = True)

      hb_tot_C      = rdfHb.Filter(selec_massfit.hb).Count().GetValue()
      bs_in_hb_C    = 0.0  
      b0_in_hb_C    = 0.0  
      bplus_in_hb_C = 0.0  
  
      if hb_tot_C > 0:
        bs_in_hb_C      = rdfHb.Filter(selec_massfit   .bs)                  .Count().GetValue() / hb_tot_C
        b0_in_hb_C      = rdfHb.Filter(selec_massfit   .b0)                  .Count().GetValue() / hb_tot_C
        bplus_in_hb_C   = rdfHb.Filter(selec_massfit   .bplus)               .Count().GetValue() / hb_tot_C

      selec_M_Hb_new_bs    = selec_M_Bs      ["phiPi_m"].Clone()
      selec_M_Hb_new_b0    = selec_M_B0      ["phiPi_m"].Clone()
      selec_M_Hb_new_bplus = selec_M_Bplus   ["phiPi_m"].Clone()

      selec_M_Hb_new_bs.Scale   ((bs_in_hb_C    / hb_tot_C)* 1.0 / selec_M_Hb_new_bs.Integral()) 
      selec_M_Hb_new_b0.Scale   ((b0_in_hb_C    / hb_tot_C)* 1.0 / selec_M_Hb_new_b0.Integral()) 
      selec_M_Hb_new_bplus.Scale((bplus_in_hb_C / hb_tot_C)* 1.0 / selec_M_Hb_new_bplus.Integral()) 
      
      selec_M_Hb_new_bs.Add(selec_M_Hb_new_b0)
      selec_M_Hb_new_bs.Add(selec_M_Hb_new_bplus)
      selec_M_Hb_new["phiPi_m"] = selec_M_Hb_new_bs

    # Data 
 
    #selec_M_Data_sf_kk   = createHistos(  baseline + kk_wrong  ,    rdfData       , variables = ["phiPi_m","dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False , sf_weights = sf_weights)
    ##selec_M_Data_sf_kk   = createHistos(  baseline + pimu_wrong      ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = True)
    #selec_M_Data_sf_pimu = createHistos(  baseline + pimu_wrong ,    rdfData       , variables = ["phiPi_m","dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, sf_weights = sf_weights)
    #selec_M_Data_sf_both = createHistos(  baseline + both_wrong ,    rdfData       , variables = ["phiPi_m","dsMu_m", "q2_lhcb_alt", "q2_coll", "kk_deltaR", "dsMu_deltaR"] , gen = False, sf_weights = sf_weights)

    selec_M_Data_sf_kk   = createHistos(  baseline + data_selec + kk_wrong   ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)
    selec_M_Data_sf_pimu = createHistos(  baseline + data_selec + pimu_wrong ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)
    #selec_M_Data_sf_both = createHistos(  baseline + data_selec + both_wrong ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)
    selec_M_Data_sf_one  = createHistos(  baseline + data_selec + wrong_sign ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)
    #selec_M_Data_sf_one  = createHistos(  baseline + data_selec + wrong_sign ,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)
    #selec_M_Data_sf_one  = createHistos(  baseline + data_selec + anti_signalRegion,    rdfData       , variables = ["phiPi_m"] , gen = False, sf_weights = sf_weights , data = True, massfit = True)

    #print(f" I have kk flip events:     {selec_M_Data_sf_kk  ['phiPi_m'].GetEntries()} ")
    #print(f" I have muflip flip events: {selec_M_Data_sf_pimu['phiPi_m'].GetEntries()} ")
    #print(f" I have both flip events:   {selec_M_Data_sf_both['phiPi_m'].GetEntries()} ")


 
    # get the signflip scale by fitting the ds mass peak of sf data against hb + signal (called hRest)

    if newHb:

      hRest       = prepareSignFlip(  selec_M_Hb_new          ["phiPi_m"].Clone()  , hb_scale_massfit,
                                     [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_M_DsTau           ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarTau       ["phiPi_m"].Clone()] ) 
   
      hRest_blind = prepareSignFlip(  selec_M_Hb_new          ["phiPi_m"].Clone()  , hb_scale_massfit,
                                     [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_M_DsTau_blind     ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarTau_blind ["phiPi_m"].Clone()] ) 
    else:

      hRest       = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                     [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_M_DsTau           ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarTau       ["phiPi_m"].Clone()] ) 
 

      hRest_blind = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                     [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_M_DsTau_blind     ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarTau_blind ["phiPi_m"].Clone()] ) 
  
    # now get the yiel ratios (postfit/prefit) 
    global scale_kk,       scale_pimu,       scale_n
    global scale_kk_blind, scale_pimu_blind, scale_n_blind
 
    kk_prefit   = selec_M_Data_sf_kk  ["phiPi_m"].Clone().Integral()
    pimu_prefit = selec_M_Data_sf_pimu["phiPi_m"].Clone().Integral() 
    bkg_prefit  = selec_M_Data_sf_one["phiPi_m"].Clone().Integral() 
    rest_prefit = hRest                          .Clone().Integral() 

    print("=============> Fit mass again to get signflip ratios")
   
    if args.findcut: 
      kk_postfit,       pimu_postfit,       rest_postfit       , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, findcut = args.findcut, cut = args.cut)
      kk_blind_postfit, pimu_blind_postfit, rest_blind_postfit , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest_blind,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, findcut = args.findcut, cut = args.cut)

    else:
      #kk_postfit,       pimu_postfit,       rest_postfit       , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3)
      #kk_blind_postfit, pimu_blind_postfit, rest_blind_postfit , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest_blind,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3)



      #hRest       = prepareSignFlip(  selec_M_Hb              ["q2_lhcb_alt"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["q2_lhcb_alt"].Clone()  , 
      #                                selec_M_DsStarMu        ["q2_lhcb_alt"].Clone()  , 
      #                                selec_M_DsTau           ["q2_lhcb_alt"].Clone()  , 
      #                                selec_M_DsStarTau       ["q2_lhcb_alt"].Clone()] ) 
 

      #fitAnotherVar(selec_M_Data_sf_kk["q2_lhcb_alt"].Clone(),selec_M_Data_sf_pimu["q2_lhcb_alt"].Clone(), selec_M_Data_sf_both["q2_lhcb_alt"].Clone(), hRest      ,selec_M_Data["q2_lhcb_alt"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "q2_lhcb_alt", -12, 12)



      #hRest       = prepareSignFlip(  selec_M_Hb              ["q2_coll"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["q2_coll"].Clone()  , 
      #                                selec_M_DsStarMu        ["q2_coll"].Clone()  , 
      #                                selec_M_DsTau           ["q2_coll"].Clone()  , 
      #                                selec_M_DsStarTau       ["q2_coll"].Clone()] ) 
 

      #fitAnotherVar(selec_M_Data_sf_kk["q2_coll"].Clone(),selec_M_Data_sf_pimu["q2_coll"].Clone(), selec_M_Data_sf_both["q2_coll"].Clone(), hRest      ,selec_M_Data["q2_coll"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "q2_coll", -12, 12)



      hRest       = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                     [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_M_DsTau           ["phiPi_m"].Clone()  , 
                                      selec_M_DsStarTau       ["phiPi_m"].Clone()] ) 
 

      #bkg_postfit,       rest_postfit, abc = fitAnotherVar(selec_M_Data_sf_one["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_m", 1.91,2.028)
      kk_postfit,       pimu_postfit,       rest_postfit       , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, findcut = args.findcut, cut = args.cut)


      #hRest       = prepareSignFlip(  selec_M_Hb              ["phiPi_deltaR"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["phiPi_deltaR"].Clone()  , 
      #                                selec_M_DsStarMu        ["phiPi_deltaR"].Clone()  , 
      #                                selec_M_DsTau           ["phiPi_deltaR"].Clone()  , 
      #                                selec_M_DsStarTau       ["phiPi_deltaR"].Clone()] ) 
 

      #kk_postfit,       pimu_postfit,       rest_postfit      =  fitAnotherVar(selec_M_Data_sf_kk["phiPi_deltaR"].Clone(),selec_M_Data_sf_pimu["phiPi_deltaR"].Clone(), selec_M_Data_sf_both["phiPi_deltaR"].Clone(), hRest      ,selec_M_Data["phiPi_deltaR"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_deltaR",0,1.2)





      #hRest       = prepareSignFlip(  selec_M_Hb              ["dsMu_m"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["dsMu_m"].Clone()  , 
      #                                selec_M_DsStarMu        ["dsMu_m"].Clone()  , 
      #                                selec_M_DsTau           ["dsMu_m"].Clone()  , 
      #                                selec_M_DsStarTau       ["dsMu_m"].Clone()] ) 
 

      #kk_postfit,       pimu_postfit,       rest_postfit      =  fitAnotherVar(selec_M_Data_sf_kk["dsMu_m"].Clone(),selec_M_Data_sf_pimu["dsMu_m"].Clone(), selec_M_Data_sf_both["dsMu_m"].Clone(), hRest      ,selec_M_Data["dsMu_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "dsMu_m",2,8)

      #kk_prefit   = selec_M_Data_sf_kk  ["dsMu_m"].Clone().Integral()
      #pimu_prefit = selec_M_Data_sf_pimu["dsMu_m"].Clone().Integral() 
      #rest_prefit = hRest                          .Clone().Integral() 

      #hRest       = prepareSignFlip(  selec_M_Hb              ["dsMu_deltaR"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["dsMu_deltaR"].Clone()  , 
      #                                selec_M_DsStarMu        ["dsMu_deltaR"].Clone()  , 
      #                                selec_M_DsTau           ["dsMu_deltaR"].Clone()  , 
      #                                selec_M_DsStarTau       ["dsMu_deltaR"].Clone()] ) 
 

      #fitAnotherVar(selec_M_Data_sf_kk["dsMu_deltaR"].Clone(),selec_M_Data_sf_pimu["dsMu_deltaR"].Clone(), selec_M_Data_sf_both["dsMu_deltaR"].Clone(), hRest      ,selec_M_Data["dsMu_deltaR"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "dsMu_deltaR", 0, 0.8)


      #hRest       = prepareSignFlip(  selec_M_Hb              ["kk_deltaR"].Clone()  , hb_scale_massfit,
      #                               [selec_M_DsMu            ["kk_deltaR"].Clone()  , 
      #                                selec_M_DsStarMu        ["kk_deltaR"].Clone()  , 
      #                                selec_M_DsTau           ["kk_deltaR"].Clone()  , 
      #                                selec_M_DsStarTau       ["kk_deltaR"].Clone()] ) 
 

      #fitAnotherVar(selec_M_Data_sf_kk["kk_deltaR"].Clone(),selec_M_Data_sf_pimu["kk_deltaR"].Clone(), selec_M_Data_sf_both["kk_deltaR"].Clone(), hRest      ,selec_M_Data["kk_deltaR"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "kk_deltaR", 0, 0.3)




    scale_kk         = kk_postfit         / kk_prefit 
    #scale_kk_blind   = kk_blind_postfit   / kk_prefit 

    scale_pimu       = pimu_postfit       / pimu_prefit 
    #scale_pimu_blind = pimu_blind_postfit / pimu_prefit 

    #scale_both       = both_postfit       / both_prefit 
    #scale_both_blind = both_blind_postfit / both_prefit 

    scale_rest       = rest_postfit       / rest_prefit 
    #scale_rest_blind = rest_blind_postfit / rest_prefit 

    #scale_bkg        = bkg_postfit       / bkg_prefit 
    #scale_bkg_blind  = bkg_blind_postfit / bkg_prefit 

    #print("=============> kk, pimu and rest scales from signflip fit are: ", scale_kk, scale_pimu, scale_rest)




  if controlPlotsHighMass:
    # if we want to do control plots, we want to plot in the high mass region 
    # the splitter must be dsMu mass and the regions is [bs mass, 8]      
    splitter     = "dsMu_m"
    regions      = [[0,8],[bsMass_, 8]]
    low_mass     = "" # remove since we are now in the high mass region
    signalRegion = "" # remove to keep as much background as possible 
 
  if controlPlotsSidebands:
    # if we want to do control plots, we want to plot in the sideband region 
    # the splitter must be phiPi mass and the regions leftSB and rightSB      
    splitter     = "phiPi_m"
    regions      = [[mlow3, mhigh3], [mlow3,mlow2], [mhigh2,mhigh3]]
    print("Doing control plots, change regions to:", regions)
    signalRegion = "" # remove since we are now in the sidebands 
    low_mass     = "" # remove to keep as much background as possible  


  ##############################
  # Start looping over regions #
  ##############################

  for i,region in enumerate(regions):

    # now we restrict to the splitter region we want to plot

    if isinstance(region, list): 
      #if region is of type list, i.e. [0,2], define selection as:
      baseline_region   = baseline + f" && ({splitter} > {region[0]}) && ({splitter} < {region[1]})" 
    else:
      #else it is a string, (used for 2D), f.e. region = "(q2_coll < 6) && (class == 0)"
      baseline_region   = baseline + f" && {region}" 
    


    # Define selections which hold all important selections next to baseline  
    selection        = baseline_region + right_sign + low_mass #for all other variables
    selection_H      = baseline_region + right_sign            #for the dsMu mass 

    selection      += score_cut  
    selection_H    += score_cut  
  
    selec            = selections(selection)
    selec_H          = selections(selection_H)
  
    ################################################
    # now start with splitter plots                #
    ################################################ 

    print(f"======> Analyzing {region} region of {splitter}")
   
 
    if not constrained:
  
      #massBins = models["phiPi_m"][0].fNbinsX

      print("=============> Fit A;B;C;D to do sideband method")
      #repeat for low mass region only
      A, B, C, S              = getABCS( rdfData, selec.bare,   "phiPi_m", sigma,   h,  dt,  binsFake = massBins, nSig = nSignalRegion, nSb = nSidebands, width = sbWidth)
    fakemass = np.genfromtxt('mass_bincontent.csv', delimiter=',')
  
    hb_tot        = rdfHb.Filter(selec.hb   + signalRegion).Count().GetValue()
    #hb_tot_C      = rdfHb.Filter(selec.hb                 ).Count().GetValue()
    hb_tot_H      = rdfHb.Filter(selec_H.hb               ).Count().GetValue()

    bs_in_hb      = 0.0  
    #bs_in_hb_C    = 0.0  
    bs_in_hb_H    = 0.0  
    b0_in_hb      = 0.0  
    #b0_in_hb_C    = 0.0  
    b0_in_hb_H    = 0.0  
    bplus_in_hb   = 0.0  
    #bplus_in_hb_C = 0.0  
    bplus_in_hb_H = 0.0  

    if hb_tot > 0:
      bs_in_hb        = rdfHb.Filter(selec   .bs + signalRegion)   .Count().GetValue() / hb_tot
      #bs_in_hb_C      = rdfHb.Filter(selec   .bs)                  .Count().GetValue() / hb_tot_C
      bs_in_hb_H      = rdfHb.Filter(selec_H .bs)                  .Count().GetValue() / hb_tot_H

      b0_in_hb        = rdfHb.Filter(selec   .b0 + signalRegion)   .Count().GetValue() / hb_tot
      #b0_in_hb_C      = rdfHb.Filter(selec   .b0)                  .Count().GetValue() / hb_tot_C
      b0_in_hb_H      = rdfHb.Filter(selec_H .b0)                  .Count().GetValue() / hb_tot_H

      bplus_in_hb     = rdfHb.Filter(selec   .bplus + signalRegion).Count().GetValue() / hb_tot
      #bplus_in_hb_C   = rdfHb.Filter(selec   .bplus)               .Count().GetValue() / hb_tot_C
      bplus_in_hb_H   = rdfHb.Filter(selec_H .bplus)               .Count().GetValue() / hb_tot_H


    ###################################################
    # Histograms in signal region (for most variables #
    ###################################################
  
    print("===> filling signal region ...")

    #debug
    #plt.figure()
    #test = rdfData.Filter(selec.bare + signalRegion).AsNumpy()
    #plt.hist(test["sf_weights"], bins=50, range=[0,50])
    #plt.savefig(f"sf_weights_ch{i}.pdf")
    #plt.clf()


    #hb_scale_S = getHbScale(selec.hb + signalRegion, selec.dsMu + signalRegion)

    #for all variables except mass plot only signal region (indicated with 'S')
    print("hammer_central, hammer_sys are:", hammer_central, hammer_sys)
    selec_S_DsMu           = createHistos(selec.dsMu        + signalRegion,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsmu"     , region = i)
    selec_S_DsMu_woHammer  = createHistos(selec.dsMu        + signalRegion,    rdfSig        , gen = False,                                                                               region = i) 
    selec_S_DsTau          = createHistos(selec.dsTau       + signalRegion,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dstau"    , region = i)
    selec_S_DsStarMu       = createHistos(selec.dsStarMu    + signalRegion,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsstarmu" , region = i)
    selec_S_DsStarTau      = createHistos(selec.dsStarTau   + signalRegion,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsstartau", region = i)
    selec_S_Hb             = createHistos(selec.hb          + signalRegion,    rdfHb         , gen = False,                                                                               region = i)
    selec_S_Mu_in_Hb       = createHistos(selec.dsMu        + signalRegion,    rdfHb         , gen = False,                                                                               region = i)
    selec_S_Data           = createHistos(selec.bare        + signalRegion,    rdfData       , gen = False,                                                                               region = i, data = True)

    #selec_S_DsTau["q2_coll"].Scale(0.3)
    #selec_S_DsStarTau["q2_coll"].Scale(0.3)


    print("I have: ", selec_S_Data["phiPi_m"].Integral(), " data events")

    # scale the tau histos to get a blind option for data fits
    selec_S_DsTau_blind     = { key: selec_S_DsTau[key].Clone()                   for key in selec_S_DsTau.keys()     }
    selec_S_DsStarTau_blind = { key: selec_S_DsStarTau[key].Clone()               for key in selec_S_DsStarTau.keys() }

    hb_scale_S = selec_S_DsMu_woHammer["phiPi_m"].Integral() * hb_ratio_massfit



    for key in selec_S_DsTau_blind.keys()    : selec_S_DsTau_blind[key]    .Scale(blind_scalar)    
    for key in selec_S_DsStarTau_blind.keys(): selec_S_DsStarTau_blind[key].Scale(blind_vector)    
 
    if newHb:
      selec_S_Bs           = createHistos(selec.bs          + signalRegion,    rdfBs         , gen = False, region = i)
      selec_S_B0           = createHistos(selec.b0          + signalRegion,    rdfB0         , gen = False, region = i)
      selec_S_Bplus        = createHistos(selec.bplus       + signalRegion,    rdfBplus      , gen = False, region = i)
  
    if constrained: 

      selec_S_Data_sf_kk   = createHistos(baseline_region + data_selec + score_cut + kk_wrong   + low_mass + signalRegion,  rdfData       , gen = False, sf_weights = sf_weights, region = i , data = True)
      selec_S_Data_sf_pimu = createHistos(baseline_region + data_selec + score_cut + pimu_wrong + low_mass + signalRegion,  rdfData       , gen = False, sf_weights = sf_weights, region = i , data = True)
      #selec_S_Data_sf_both = createHistos(baseline_region + data_selec + score_cut + both_wrong + low_mass + signalRegion,  rdfData       , gen = False, sf_weights = sf_weights, region = i , data = True)
      #selec_S_Data_sf_one  = createHistos(baseline_region + data_selec + score_cut + wrong_sign + low_mass + signalRegion,  rdfData       , gen = False, sf_weights = sf_weights, region = i , data = True)

      #selec_S_DataL      = createHistos(selec.bare        + leftSB  + score_cut,         rdfData       , gen = False, region = i)
      #selec_S_DataR      = createHistos(selec.bare        + rightSB + score_cut,         rdfData       , gen = False, region = i)
   
    print("===> signal region done...")
  
    print("===> filling total region ...")
    ###################################################
    # Histograms in complete region (for Ds mass)     #
    ###################################################
 
    """
    hb_scale_C = getHbScale(selec.hb + signalRegion, selec.dsMu + signalRegion)
  
    #for all variables except mass plot only signal region (indicated with 'S')
    selec_C_DsMu           = createHistos(selec.dsMu        ,    rdfSig        , gen = False, histSys = addSys, sig = "dsmu"     )
    selec_C_DsTau          = createHistos(selec.dsTau       ,    rdfSig        , gen = False, histSys = addSys, sig = "dstau"    )
    selec_C_DsStarMu       = createHistos(selec.dsStarMu    ,    rdfSig        , gen = False, histSys = addSys, sig = "dsstarmu" )
    selec_C_DsStarTau      = createHistos(selec.dsStarTau   ,    rdfSig        , gen = False, histSys = addSys, sig = "dsstartau")
    print("=> signals done ...")
    selec_C_Hb             = createHistos(selec.hb          ,    rdfHb         , gen = False)
    selec_C_Data           = createHistos(selec.bare        ,    rdfData       , gen = False)

    # scale the tau histos to get a blind option for data fits
    selec_C_DsTau_blind     = { key: selec_C_DsTau[key].Clone()                   for key in selec_C_DsTau.keys()     }
    selec_C_DsStarTau_blind = { key: selec_C_DsStarTau[key].Clone()               for key in selec_C_DsStarTau.keys() }

    for key in selec_C_DsTau_blind.keys()    : selec_C_DsTau_blind[key]    .Scale(blind_scalar)    
    for key in selec_C_DsStarTau_blind.keys(): selec_C_DsStarTau_blind[key].Scale(blind_vector)    
 
    if newHb:
      selec_C_Bs           = createHistos(selec.bs          ,    rdfBs         , gen = False)
      selec_C_B0           = createHistos(selec.b0          ,    rdfB0         , gen = False)
      selec_C_Bplus        = createHistos(selec.bplus       ,    rdfBplus      , gen = False)
  
    if constrained: 

      #print("using baseline; ", baseline_region)
      selec_C_Data_sf_kk   = createHistos(  baseline_region + score_cut + kk_wrong   + low_mass ,  rdfData       , gen = False)
      selec_C_Data_sf_pimu = createHistos(  baseline_region + score_cut + pimu_wrong + low_mass ,  rdfData       , gen = False)
  
      #get the signflip scale by fitting the ds mass peak of sf data against hb + signal (called hRest)
      hRest       = prepareSignFlip(  selec_C_Hb              ["phiPi_m"].Clone()  , 
                                      hb_scale_C                                   ,
                                     [selec_C_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_C_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_C_DsTau           ["phiPi_m"].Clone()  , 
                                      selec_C_DsStarTau       ["phiPi_m"].Clone()]  ) 
  
      hRest_blind = prepareSignFlip(  selec_C_Hb              ["phiPi_m"].Clone()  , 
                                      hb_scale_C                                   ,
                                     [selec_C_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_C_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_C_DsTau_blind     ["phiPi_m"].Clone()  , 
                                      selec_C_DsStarTau_blind ["phiPi_m"].Clone()]  ) 
 
      print("phipi mu histo has events: ",   selec_C_DsMu            ["phiPi_m"].Clone().Integral() )
      print("phipi mu* histo has events: ",  selec_C_DsStarMu            ["phiPi_m"].Clone().Integral() )
      print("phipi tau histo has events: ",  selec_C_DsTau            ["phiPi_m"].Clone().Integral() )
      print("phipi tau* histo has events: ", selec_C_DsStarTau           ["phiPi_m"].Clone().Integral() )
 
      global scale_C_kk,       scale_C_bkg,       scale_C_n
      global scale_C_kk_blind, scale_C_bkg_blind, scale_C_n_blind
  
      print("=============> Fit mass again to get signflip ratios")
      scale_C_kk,       scale_C_bkg,       scale_C_n       = getSignflipRatio(selec_C_Data_sf_kk["phiPi_m"].Clone(),selec_C_Data_sf_pimu["phiPi_m"].Clone(),hRest      ,selec_C_Data["phiPi_m"].Clone())
      scale_C_kk_blind, scale_C_bkg_blind, scale_C_n_blind = getSignflipRatio(selec_C_Data_sf_kk["phiPi_m"].Clone(),selec_C_Data_sf_pimu["phiPi_m"].Clone(),hRest_blind,selec_C_Data["phiPi_m"].Clone())

      print("sclaes from fit are: ", scale_C_kk, scale_C_bkg, scale_C_n)


    else:
        selec_C_DataL        = createHistos(selec.bare        + leftSB,          rdfData       , gen = False)
        selec_C_DataR        = createHistos(selec.bare        + rightSB,         rdfData       , gen = False)
  



    ########################################################
    # Histograms in high&low mass region (For Ds+mu mass)  #
    ########################################################
  
    hb_scale_H = getHbScale(selec.hb + signalRegion, selec.dsMu + signalRegion)
  
    #for all variables except mass plot only signal region (indicated with 'S')
    selec_H_DsMu           = createHistos(selec_H.dsMu        ,    rdfSig        , gen = False, histSys = addSys, sig = "dsmu"     )
    selec_H_DsTau          = createHistos(selec_H.dsTau       ,    rdfSig        , gen = False, histSys = addSys, sig = "dstau"    )
    selec_H_DsStarMu       = createHistos(selec_H.dsStarMu    ,    rdfSig        , gen = False, histSys = addSys, sig = "dsstarmu" )
    selec_H_DsStarTau      = createHistos(selec_H.dsStarTau   ,    rdfSig        , gen = False, histSys = addSys, sig = "dsstartau")
    print("=> signals done ...")
    selec_H_Hb             = createHistos(selec_H.hb          ,    rdfHb         , gen = False)
    selec_H_Data           = createHistos(selec_H.bare        ,    rdfData       , gen = False)

    # scale the tau histos to get a blind option for data fits
    selec_H_DsTau_blind     = { key: selec_H_DsTau[key].Clone()                   for key in selec_H_DsTau.keys()     }
    selec_H_DsStarTau_blind = { key: selec_H_DsStarTau[key].Clone()               for key in selec_H_DsStarTau.keys() }

    for key in selec_H_DsTau_blind.keys()    : selec_H_DsTau_blind[key]    .Scale(blind_scalar)    
    for key in selec_H_DsStarTau_blind.keys(): selec_H_DsStarTau_blind[key].Scale(blind_vector)    
 
    if newHb:
      selec_H_Bs           = createHistos(selec_H.bs          ,    rdfBs         , gen = False)
      selec_H_B0           = createHistos(selec_H.b0          ,    rdfB0         , gen = False)
      selec_H_Bplus        = createHistos(selec_H.bplus       ,    rdfBplus      , gen = False)
  
    if constrained: 

      #print("using baseline; ", baseline_region)
      selec_H_Data_sf_kk   = createHistos(  baseline_region + score_cut + kk_wrong   ,  rdfData       , gen = False)
      selec_H_Data_sf_pimu = createHistos(  baseline_region + score_cut + pimu_wrong ,  rdfData       , gen = False)
  
      #get the signflip scale by fitting the ds mass peak of sf data against hb + signal (called hRest)
      hRest       = prepareSignFlip(  selec_H_Hb              ["phiPi_m"].Clone()  , 
                                      hb_scale_H                                   ,
                                     [selec_H_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_H_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_H_DsTau           ["phiPi_m"].Clone()  , 
                                      selec_H_DsStarTau       ["phiPi_m"].Clone()]  ) 
  
      hRest_blind = prepareSignFlip(  selec_H_Hb              ["phiPi_m"].Clone()  , 
                                      hb_scale_H                                   ,
                                     [selec_H_DsMu            ["phiPi_m"].Clone()  , 
                                      selec_H_DsStarMu        ["phiPi_m"].Clone()  , 
                                      selec_H_DsTau_blind     ["phiPi_m"].Clone()  , 
                                      selec_H_DsStarTau_blind ["phiPi_m"].Clone()]  ) 
 
      print("phipi mu histo has events: ",   selec_H_DsMu            ["phiPi_m"].Clone().Integral() )
      print("phipi mu* histo has events: ",  selec_H_DsStarMu            ["phiPi_m"].Clone().Integral() )
      print("phipi tau histo has events: ",  selec_H_DsTau            ["phiPi_m"].Clone().Integral() )
      print("phipi tau* histo has events: ", selec_H_DsStarTau           ["phiPi_m"].Clone().Integral() )
 
      global scale_H_kk,       scale_H_bkg,       scale_H_n
      global scale_H_kk_blind, scale_H_bkg_blind, scale_H_n_blind
  
      print("=============> Fit mass again to get signflip ratios")
      scale_H_kk,       scale_H_bkg,       scale_H_n       = getSignflipRatio(selec_H_Data_sf_kk["phiPi_m"].Clone(),selec_H_Data_sf_pimu["phiPi_m"].Clone(),hRest      ,selec_H_Data["phiPi_m"].Clone())
      scale_H_kk_blind, scale_H_bkg_blind, scale_H_n_blind = getSignflipRatio(selec_H_Data_sf_kk["phiPi_m"].Clone(),selec_H_Data_sf_pimu["phiPi_m"].Clone(),hRest_blind,selec_H_Data["phiPi_m"].Clone())

      print("sclaes from fit are: ", scale_H_kk, scale_H_bkg, scale_H_n)


    else:
        selec_H_DataL        = createHistos(selec_H.bare        + leftSB,          rdfData       , gen = False)
        selec_H_DataR        = createHistos(selec_H.bare        + rightSB,         rdfData       , gen = False)
  
    """
    
    print("===> high mass region done...")

    if constrained: name = "constrained"
    else: name = "unconstrained"
  
    if newHb: name += "_newHb"
    else: name += ""
  
    if pastNN: name += "_pastNN"
    else: name += ""
  
    for var in models.keys():
   
      #if (var == splitter): continue; 
 
      #create root file for every variable which holds shapes
      shapes_name = f"/{var}_shapes_{name}_{splitter}_ch{i}.root"
      
      myFile       = ROOT.TFile.Open(shapes_folder       + shapes_name, "RECREATE")
      myFile_blind = ROOT.TFile.Open(shapes_folder_blind + shapes_name, "RECREATE")
      print(f"===> Producing stacked plot for variable: {var}") 
      if "gen" in var:
        #skip gen variables
        print("This is a gen variable... skip!")
        continue
   
      ######################################
      ## Plot variables except phiPi mass ##
      ######################################
  
      if  (var!= "" ):
   


 
        histos = {"dsTau"    :  selec_S_DsTau[var]    .Clone(),
                 "dsStarTau" :  selec_S_DsStarTau[var].Clone(), 
                 "dsMu"      :  selec_S_DsMu[var]     .Clone(), 
                 "dsStarMu"  :  selec_S_DsStarMu[var] .Clone(), 
                 "hb"        :  selec_S_Hb[var]       .Clone(), 
                 "data"      :  selec_S_Data[var]     .Clone()}
  
  
        if newHb:
          histos["bs"]    = selec_S_Bs[var]           .Clone() 
          histos["b0"]    = selec_S_B0[var]           .Clone() 
          histos["bplus"] = selec_S_Bplus[var]        .Clone() 
  

        # Remark: histos_blind = histos.copy() only does a shallow copy! It defines a new
        # address for the dict, but not for the elements stored in the dict. Save way is:

        histos_blind = {key: hist.Clone() for key, hist in histos.items()}
        histos_blind["dsTau"]     = selec_S_DsTau_blind    [var].Clone()   
        histos_blind["dsStarTau"] = selec_S_DsStarTau_blind[var].Clone()   
 
        print("first check here, data blind: ", histos_blind["data"].Integral() , "data unblinded:", histos["data"].Integral() )
 
        if hammer_sys: 
          histos       = addSystematics(histos      , var, selec_S_DsTau      , selec_S_DsMu, selec_S_DsStarTau      , selec_S_DsStarMu )
          histos_blind = addSystematics(histos_blind, var, selec_S_DsTau_blind, selec_S_DsMu, selec_S_DsStarTau_blind, selec_S_DsStarMu )

        if constrained:
 
          histos      ["data_sf_kk"]   = selec_S_Data_sf_kk[var]  .Clone() 
          histos      ["data_sf_pimu"] = selec_S_Data_sf_pimu[var].Clone() 
          #histos      ["data_sf_both"] = selec_S_Data_sf_both[var].Clone() 
          #histos      ["data_sf_one"]  = selec_S_Data_sf_one[var].Clone() 
          #histos_blind["data_sf_kk"]   = selec_S_Data_sf_kk[var]  .Clone() 
          #histos_blind["data_sf_pimu"] = selec_S_Data_sf_pimu[var].Clone() 
          #histos_blind["data_sf_both"] = selec_S_Data_sf_both[var].Clone() 

          #histos      ["combL"] = selec_S_DataL[var].Clone() 
          #histos      ["combR"] = selec_S_DataR[var].Clone() 
          #histos_blind["combL"] = selec_S_DataL[var].Clone() 
          #histos_blind["combR"] = selec_S_DataR[var].Clone() 

          #A = abc[0]
          #B = abc[1]
          #C = abc[2]

          #histosScaled       , _ = stackedPlot(histos,       var, hb_scale_S, fakemass = fakemass, abc = abc, scale_bkg = scale_bkg, scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          histosScaled       , _ = stackedPlot(histos,       var, hb_scale_S, abc = abc,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          #histosScaledLog    , _ = stackedPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False , log = True)
          #histosScaled_blind , _ = stackedPlot(histos_blind, var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk_blind, scale_pimu = scale_pimu_blind, scale_rest = scale_rest_blind, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True  )
 
          #histosScaled_shapes    = shapesPlot(histos,        var, hb_scale_S, fakemass = fakemass, abc = abc, scale_bkg = scale_bkg, scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          #histosScaled_shapes    = shapesPlot(histos,        var, hb_scale_S, fakemass = fakemass, abc = abc, scale_bkg = scale_bkg, scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          #histosScaled_shapes = shapesPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          #histosScaled_shapes_log = shapesPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False , log = True)
       

          #histosScaled_shapes_bkg = shapesBkg(histos,       var, region = f"{splitter}_ch{i}")
 
 
        else:
 
          histos      ["combL"] = selec_S_DataL[var].Clone() 
          histos      ["combR"] = selec_S_DataR[var].Clone() 
          histos_blind["combL"] = selec_S_DataL[var].Clone() 
          histos_blind["combR"] = selec_S_DataR[var].Clone() 
  
          histosScaled       , _ = stackedPlot(histos,       var, hb_scale_S, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False)
          histosScaled_blind , _ = stackedPlot(histos_blind, var, hb_scale_S, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True )
  
 
        writeShapes(histosScaled,       myFile,       binned = True, channel = i)
        #writeShapes(histosScaled_blind, myFile_blind, binned = True, channel = i)

        card       = writeBinnedDatacard(histosScaled      , var, i, splitter              )
        #card_blind = writeBinnedDatacard(histosScaled_blind, var, i, splitter, blind = True)

        cards[var].append(card)
        #cards_blind[var].append(card_blind)


        #find optimal working points
        if "WP" in models[var][0].fName:
          dstau   = [] 
          dssttau = [] 
          dsmu    = [] 
          dsstmu  = [] 
          hb      = [] 
          comb    = [] 
          #take one to get the binning
          bins  = histosScaled["dsMu"].GetNbinsX()  #f.e. 100
 
          #reproduce this binning with numpy
          lower = models[var][0].fXLow #0.0
          upper = models[var][0].fXUp  #f.e. 0.1
          #attention, linspace expects the nr of steps, not the bins!
          #steps = bins + 1!
          edges = np.linspace(lower,upper,bins+1)#f.e. [0,0.001,0.002,...]

          #this is a finely binned histogram over the scores, lets count S and B in every bin!
          for i in range(bins):
            dstau   .append(histosScaled["dsTau"]    .GetBinContent(i+1))
            dssttau .append(histosScaled["dsStarTau"].GetBinContent(i+1))
            #dsmu    .append(histosScaled["dsMu"]     .GetBinContent(i+1))
            #dsstmu  .append(histosScaled["dsStarMu"] .GetBinContent(i+1))
            #hb      .append(histosScaled["hb"]       .GetBinContent(i+1))
            comb    .append(histosScaled["comb"]     .GetBinContent(i+1))

          S = []
          B = []

          if "lower" in models[var][0].fName:
            #this is a cut where we keep everything "lower" than the cut
            for i in range(bins):
            
              # up until bin i ( in nclusive i) we have signals:
              # example: bin content is b = [1,2,3,4], then: b[:1] = [1], b[:2] = [1,2], etc.
              S.append(sum(dstau[:i+1]) + sum(dssttau[:i+1])) #+1 since i starts at 0                     
              B.append(sum(comb[:i+1]))  

            cut = edges[1:] 

          if "upper" in models[var][0].fName:
            #this is a cut where we keep everything "higher" than the cut
            for i in range(bins):
            
              S.append(sum(dstau[i:]) + sum(dssttau[i:])) #now we start at 0                     
              B.append(sum(comb[i:]))  


            cut = edges[: len(edges)-1]

          fom = [s / np.sqrt(b) for (s,b) in zip(S,B)]

          #find cut at which fom is max
          max_fom = np.max(fom)
          max_cut = cut[np.argmax(fom)]


          plt.figure()
          plt.plot(cut,fom)
          plt.scatter(cut,fom)
          plt.xlabel(f"Cut on {var}") 
          plt.ylabel(r"$ S / \sqrt{B}$ ")
          plt.ylim(0, max_fom * 1.2)
          plt.text(max_cut, max_fom * 1.1, f"Optimal cut at {max_cut}")
          toSave = f"/work/pahwagne/RDsTools/plots/cmsplots_binned/{dt}/"
          plt.savefig(toSave + f"{var}_WP.pdf")
          plt.savefig(toSave + f"{var}_WP.png")
          plt.clf()


      ######################################
      ## phiPi mass ##
      ######################################
      """
      if  (var == "phiPi_m"):
    
 
        histos = {"dsTau"    :  selec_C_DsTau[var]    .Clone(),
                 "dsStarTau" :  selec_C_DsStarTau[var].Clone(), 
                 "dsMu"      :  selec_C_DsMu[var]     .Clone(), 
                 "dsStarMu"  :  selec_C_DsStarMu[var] .Clone(), 
                 "hb"        :  selec_C_Hb[var]       .Clone(), 
                 "data"      :  selec_C_Data[var]     .Clone()}
  
  
        if newHb:
          histos["bs"]    = selec_C_Bs[var]           .Clone() 
          histos["b0"]    = selec_C_B0[var]           .Clone() 
          histos["bplus"] = selec_C_Bplus[var]        .Clone() 
  

        # Remark: histos_blind = histos.copy() only does a shallow copy! It defines a new
        # address for the dict, but not for the elements stored in the dict. Save way is:

        histos_blind = {key: hist.Clone() for key, hist in histos.items()}
        histos_blind["dsTau"]     = selec_C_DsTau_blind    [var].Clone()   
        histos_blind["dsStarTau"] = selec_C_DsStarTau_blind[var].Clone()   
 
        print("first check here, data blind: ", histos_blind["data"].Integral() , "data unblinded:", histos["data"].Integral() )
 
        if addSys: 
          histos       = addSystematics(histos      , var, selec_C_DsTau, selec_C_DsMu, selec_C_DsStarTau, selec_C_DsStarMu )
          histos_blind = addSystematics(histos_blind, var, selec_C_DsTau, selec_C_DsMu, selec_C_DsStarTau, selec_C_DsStarMu )

        if constrained:
 
          histos      ["data_sf_kk"]   = selec_C_Data_sf_kk[var]  .Clone() 
          histos      ["data_sf_pimu"] = selec_C_Data_sf_pimu[var].Clone() 
          histos_blind["data_sf_kk"]   = selec_C_Data_sf_kk[var]  .Clone() 
          histos_blind["data_sf_pimu"] = selec_C_Data_sf_pimu[var].Clone() 

          #histos      ["combL"] = selec_C_DataL[var].Clone() 
          #histos      ["combR"] = selec_C_DataR[var].Clone() 
          #histos_blind["combL"] = selec_C_DataL[var].Clone() 
          #histos_blind["combR"] = selec_C_DataR[var].Clone() 


          A = abc[0]
          B = abc[1]
          C = abc[2]
          histosScaled       , _ = stackedPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          histosScaled_blind , _ = stackedPlot(histos_blind, var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk_blind, scale_pimu = scale_pimu_blind, scale_rest = scale_rest_blind, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True  )
 
          histosScaled_shapes = shapesPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          histosScaled_shapes_log = shapesPlot(histos,       var, hb_scale_S, A = A, B = B, C = C,  scale_kk = scale_kk,       scale_pimu = scale_pimu,       scale_rest = scale_rest, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False , log = True)
        
 
        else:
 
          histos      ["combL"] = selec_C_DataL[var].Clone() 
          histos      ["combR"] = selec_C_DataR[var].Clone() 
          histos_blind["combL"] = selec_C_DataL[var].Clone() 
          histos_blind["combR"] = selec_C_DataR[var].Clone() 
  
          histosScaled       , _ = stackedPlot(histos,       var, hb_scale_S, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False)
          histosScaled_blind , _ = stackedPlot(histos_blind, var, hb_scale_S, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True )
  
 
        writeShapes(histosScaled,       myFile,       binned = True, channel = i)
        writeShapes(histosScaled_blind, myFile_blind, binned = True, channel = i)

        card       = writeBinnedDatacard(histosScaled      , var, i, splitter              )
        card_blind = writeBinnedDatacard(histosScaled_blind, var, i, splitter, blind = True)

        cards[var].append(card)
        cards_blind[var].append(card_blind)

      """
      ######################################
      ## Dsmu mass ##
      ######################################

      """
      if (var == "dsMu_m") :
    
 
        histos = {"dsTau"    :  selec_H_DsTau[var]    .Clone(),
                 "dsStarTau" :  selec_H_DsStarTau[var].Clone(), 
                 "dsMu"      :  selec_H_DsMu[var]     .Clone(), 
                 "dsStarMu"  :  selec_H_DsStarMu[var] .Clone(), 
                 "hb"        :  selec_H_Hb[var]       .Clone(), 
                 "data"      :  selec_H_Data[var]     .Clone()}
  
        if newHb:
          histos["bs"]    = selec_H_Bs[var]           .Clone() 
          histos["b0"]    = selec_H_B0[var]           .Clone() 
          histos["bplus"] = selec_H_Bplus[var]        .Clone() 
  

        # Remark: histos_blind = histos.copy() only does a shallow copy! It defines a new
        # address for the dict, but not for the elements stored in the dict. Save way is:

        histos_blind = {key: hist.Clone() for key, hist in histos.items()}
        histos_blind["dsTau"]     = selec_H_DsTau_blind    [var].Clone()   
        histos_blind["dsStarTau"] = selec_H_DsStarTau_blind[var].Clone()   
 
        print("first check here, data blind: ", histos_blind["data"].Integral() , "data unblinded:", histos["data"].Integral() )
 
        if addSys: 
          histos       = addSystematics(histos      , var, selec_H_DsTau, selec_H_DsMu, selec_H_DsStarTau, selec_H_DsStarMu )
          histos_blind = addSystematics(histos_blind, var, selec_H_DsTau, selec_H_DsMu, selec_H_DsStarTau, selec_H_DsStarMu )

        if constrained:
 
          histos      ["data_sf_kk"]   = selec_H_Data_sf_kk[var]  .Clone() 
          histos      ["data_sf_pimu"] = selec_H_Data_sf_pimu[var].Clone() 
          histos_blind["data_sf_kk"]   = selec_H_Data_sf_kk[var]  .Clone() 
          histos_blind["data_sf_pimu"] = selec_H_Data_sf_pimu[var].Clone() 

          print("Before calling stackedPlot: Number of tau* signals: ", histos["dsStarTau"].Integral(), " Number of mu* signals: ", histos["dsStarMu"].Integral())
          print("Before calling stackedPlot: Blind number of tau* signals: ", histos_blind["dsStarTau"].Integral(), " Number of mu* signals: ", histos_blind["dsStarMu"].Integral())

          print("then check here, data blind: ", histos_blind["data"].Integral() , "data unblinded:", histos["data"].Integral() )
 
          histosScaled       = stackedPlot(histos,       var, hb_scale_H, scale_kk = scale_H_kk,       scale_bkg = scale_H_bkg,       scale_n = scale_H_n,       bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False )
          histosScaled_blind = stackedPlot(histos_blind, var, hb_scale_H, scale_kk = scale_H_kk_blind, scale_bkg = scale_H_bkg_blind, scale_n = scale_H_n_blind, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True  )
 
        
 
        else:
 
          histos      ["combL"] = selec_H_DataL[var].Clone() 
          histos      ["combR"] = selec_H_DataR[var].Clone() 
          histos_blind["combL"] = selec_H_DataL[var].Clone() 
          histos_blind["combR"] = selec_H_DataR[var].Clone() 
  
          histosScaled       = stackedPlot(histos,       var, hb_scale_H, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = False)
          histosScaled_blind = stackedPlot(histos_blind, var, hb_scale_H, fakemass = fakemass, A = A, B = B, C = C, S = S, bs = bs_in_hb, b0 = b0_in_hb, bplus = bplus_in_hb, region = f"{splitter}_ch{i}", blind = True )
  
  
        writeShapes(histosScaled,       myFile,       binned = True, channel = i)
        writeShapes(histosScaled_blind, myFile_blind, binned = True, channel = i)

        card       = writeBinnedDatacard(histosScaled      , var, i, splitter              )
        card_blind = writeBinnedDatacard(histosScaled_blind, var, i, splitter, blind = True)

        cards[var].append(card)
        cards_blind[var].append(card_blind)
        """




#more advanced, create binned datacards
if control == "highmass":
  createBinnedPlots(split,binning, controlPlotsHighMass = True) 
elif control == "sb":
  createBinnedPlots(split,binning, controlPlotsSidebands = True) 
else:
  createBinnedPlots(split,binning) 


