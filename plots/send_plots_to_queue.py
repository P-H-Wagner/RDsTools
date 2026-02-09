import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
import json

#ROOT.ROOT.EnableImplicitMT() 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio, getSignflipRatioTest, fitAnotherVar
from helper import * 
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels, special_models, special_models_q2_coll
from blinding import *
import numpy as np


########################
# Argument parsing     #
########################

def boolean_string(s):
    if s not in {"false", "true"}:
        raise ValueError("Error: Not a valid boolean string, please use 'true' or 'false' (all lowercase!)")
    return s == "true"

# parsing
parser = argparse.ArgumentParser()
parser.add_argument("--nn",          required = True,     help = "Specify 'before' or 'past' neural network plots ") 
parser.add_argument("--hammer",      required = True,     help = "Specify 'true' or 'false' to apply hammer weights") 
parser.add_argument("--hammer_sys",  required = True,     help = "Specify 'true' or 'false' to save weight variation shapes") 
parser.add_argument("--trigger",     required = True,     help = "Specify mu7 or mu9 for the trigger") 
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
parser.add_argument("--control",                          help = "If given, run control plots, either 'highmass', 'leftsb', 'rightsb' or 'complete' or 'custom' ") 
parser.add_argument("--cut",                              help = "Cut on the discriminator score 5. ") 
#parser.add_argument("--findcut",     action='store_true', help = "If given, we run thecut scan") 
args = parser.parse_args()

if args.nn not in  ["before", "past"]:
  raise ValueError ("Error: Not a valid key for --nn, please use 'before' or 'past' (all lowercase!)")
else: pastNN      = (args.nn == "past")

if args.hammer not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --hammer, please use 'true' or 'false' (all lowercase!)")
else: hammer_central = (args.hammer == "true")

if args.hammer_sys not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --sys, please use 'true' or 'false' (all lowercase!)")
else: hammer_sys = (args.hammer_sys == "true")

if args.trigger not in ["mu7", "mu9"]:
  raise ValueError ("Error: Not a valid key for --trigger, please use 'mu7' or 'mu9'")
else: trigger = args.trigger

if args.bdt not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt , please use 'true' or 'false' (all lowercase!)")
else: sf_weights = (args.bdt == "true")

if args.control and args.control not in ["highmass","leftsb","rightsb", "complete", "custom"]:
    raise ValueError ("Error: Not a valid key for --control, please use 'highmass', 'sb' or 'complete' or 'custom' ")
control = args.control

if args.cut  : score_cut = f" && (score5 <= {args.cut}) "
else         : score_cut = ""

if args.debug: debug = 50000
else         : debug = None

print(f"====> Running  {args.nn} neural network on trigger {trigger} and sf_weights: {sf_weights} with hammer on: {hammer_central} and sys {hammer_sys}")
if pastNN: models.update(pastNN_models)

if hammer_sys:
  print(f"====> Adding the following systematics: {systematics_scalar} and {systematics_vector}") #defined in helper.py
  sys_scalar = systematics_scalar
  sys_vector = systematics_vector

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/20_10_2025_09_57_12/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)


# Data, selections, trigger
baseline      = base_wout_tv_25
print("Using baseline selection:", baseline)

#before NN
bdt_data   = bdt_data_25
sig_cons   = sig_cons_25
hb_cons    = hb_cons_25
data_cons  = data_cons_25

#after NN
if trigger == "mu7": 
  cons_pastNN_25 = cons_pastNN_25_mu7 
else:                
  cons_pastNN_25 = cons_pastNN_25_mu9

sig_cons_pastNN     =  cons_pastNN_25["sig"]    
hb_cons_pastNN      =  cons_pastNN_25["hb"]     
data_cons_pastNN    =  cons_pastNN_25["data"]   

skimmed = True 

# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

# get date and time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S") 

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

with open( toSave_plots + f"/info.txt", "a") as f:
  f.write( f"====> Running  {args.nn} neural network on trigger {trigger} and sf_weights: {sf_weights} with hammer on: {hammer_central} and sys {hammer_sys} \n")

#############
# Split MC  #
#############

if trigger == "mu7":
  data_selec   = " && (mu7_ip4 == 1)"
  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
else:
  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "

#also freeze sigma here!
sigma = 0.009


# set the to be splitter variable and binning
#split   = "cosMuW_lhcb_alt"
split   = "q2_coll"
#split   = "q2_lhcb_alt"
#split   = "class"
#split   = "phiPi_m"
#split   = "kk_m"
#split   = "e_star_reco_weighted"
#split = "score1"
#split = "score2"
#binning = [[0,0.03],[0.03,0.06],[0.06,0.1],[0.1,0.2],[0.2,0.4],[0.4,0.6],[0.6,1.0]] #score 3
#binning = [[0,0.02],[0.02,0.1],[0.1,0.4],[0.4,1]] #score 2
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,12]]
binning = [[0,6],[6,7],[7,8],[8,9],[9,10],[10,12]]
#binning = [[0,6],[6,7],[7,8],[8,9],[9,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[8,8.5],[8.5,9],[9,10],[10,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,12]]
#binning = [[0,4],[4,5],[5,6],[6,7],[7,8],[9,12]]
#binning = [[-1,-0.5],[-0.5, -0.0], [0,0.5],[0.5,1]] #cosmuw
#binning = [[1.94, 1.95],[1.95, 1.96], [1.96, 1.98], [1.98,2.0]] #mass]
#binning = [[0,0.5], [0.5,1.0],[1.0, 1.5], [1.5, 2.0],[2.0, 3.0]] # estar
#binning = [[1.0,1.015], [1.015,1.02],[1.02, 1.025], [1.025, 1.035]] # estar
#binning = [[0,7],[7,12]]
#binning = [[-100000,1000000]]
#binning = [[-0.01,0.01],[0.99,1.01],[1.99,2.01],[2.99,3.01],[3.99,4.01],[4.99,5.01]]
#binning = [[0.99,1.01]]
#binning = [[4.99,5.01]]
binning = [[-99,99]]

################
## 3D binning ##
################

binning_str = []
#
sr_string = " && ((phiPi_m >= 1.94134) && (phiPi_m <= 1.99534 ) )"
sb_string = " && ((phiPi_m < 1.94 ) || (phiPi_m > 1.995  ) )"
#

#custom strings
mass_r1 = f" (  (phiPi_m <  {dsMass_ - 5*sigma}) || ( phiPi_m > {dsMass_ + 5*sigma} ))"  
mass_r2 = f" ( ((phiPi_m >= {dsMass_ - 5*sigma}) && ( phiPi_m < {dsMass_ - 3*sigma} )) || ((phiPi_m <= {dsMass_ + 5*sigma}) && ( phiPi_m > {dsMass_ + 3*sigma} )) )" 
mass_r3 = f" ( ((phiPi_m >= {dsMass_ - 3*sigma}) && ( phiPi_m < {dsMass_ - 2*sigma} )) || ((phiPi_m <= {dsMass_ + 3*sigma}) && ( phiPi_m > {dsMass_ + 2*sigma} )) )" 
mass_r4 = f" ( ((phiPi_m >= {dsMass_ - 2*sigma}) && ( phiPi_m < {dsMass_ - 1*sigma} )) || ((phiPi_m <= {dsMass_ + 2*sigma}) && ( phiPi_m > {dsMass_ + 1*sigma} )) )" 
mass_r5 = f" ( ((phiPi_m >= {dsMass_ - 1*sigma}) && ( phiPi_m < {dsMass_ - 0*sigma} )) || ((phiPi_m <= {dsMass_ + 1*sigma}) && ( phiPi_m > {dsMass_ + 0*sigma} )) )" 
#
mass_binning = [mass_r1, mass_r2, mass_r3, mass_r4, mass_r5]

for i,mass_b in enumerate(mass_binning):

  if i in [0,1,2,3]:
    #first mass bin, no slice in q2
    binning_str.append(mass_b)

  else: 
    for b in binning:
      binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} ) && " + mass_b)

#for b in binning:
#  
#  # b is a list
#  binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" + sr_string)  
#  #binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" + sb_string)  
#  #binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" )  
#
  

##overwrite
binning = binning_str

#print("====> Using 3D regions: ")
#print(binning_str)

##############################
# Load chain into RDataFrame #
##############################

def getRdf(dateTimes, debug = None, skimmed = None, rdfSys = False, sf_weights = None, bph_part = None, isData = None):

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
 
      #if debug and isData:
      #  print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHh")
      #  files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*22*" #test

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


    self.hb =        selec + mc_selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1)  && (gen_same_mother == 1) "
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



  sig_sample = sig_cons_pastNN
  print("sample is:", sig_cons_pastNN)

  chainSigSB, rdfSigSB     = getRdf(sig_cons_pastNN                 , debug = debug)
  chainSig,   rdfSig       = getRdf(sig_cons_pastNN                 , debug = debug)
  chainHb,    rdfHb        = getRdf(hb_cons_pastNN                  , debug = debug)
  _,          rdfHbHammer  = getRdf(hb_cons_pastNN                  , debug = debug)
  chainData,  rdfData      = getRdf(data_cons_pastNN                , debug = debug)

  print("rdf data has events: ", rdfData.Count().GetValue() )
  print("rdf sig has events: ", rdfSig.Count().GetValue() )


else:

  print("Producing pre-NN plots")

  if hammer_central: sig_cons = sig_cons_hammer_25

  chainSigSB, rdfSigSB     = getRdf(sig_cons                              , debug = debug, skimmed = skimmed, rdfSys = hammer_central )# , skimmed = baseline_name)#, debug = 1)
  chainSig,   rdfSig       = getRdf(sig_cons                              , debug = debug, skimmed = skimmed, rdfSys = hammer_central )# , skimmed = baseline_name)#, debug = 1)
  chainHb,    rdfHb        = getRdf(hb_cons                               , debug = debug, skimmed = skimmed )# , skimmed = baseline_name)#, debug = 1)
  chainData, rdfData       = getRdf(data_cons    ,sf_weights = sf_weights , debug = debug, skimmed = skimmed, isData = True )# , skimmed = baseline_name)#, debug = 1)
  #print("---------------> rdf has events: ", rdfData.Count().GetValue() )


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

  print(" ====> We are in Region", region)


  histos = {}

  if pastNN:           models.update(pastNN_models)
  if massfit == False: models.update(modelsSR)

  for var,model in models.items(): 

      #print("var is", var, "split is", split)

      ################################
      # Special binnings             #
      ################################
      if var == "score1" and (split == "class") and region != None:
        
        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")
  
        model = special_models[var + f"_bin{region}" ]

      #if var == "score1" and (split == "q2_coll" or split == "q2_lhcb_alt") and region != None:
      #  
      #  #adapt the binning
      #  print(f"Adapt binning for {var} and region {region}")
  
      #  model = special_models_q2_coll[var + f"_bin{region}" ]                                     


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
          histos[var] = rdf.Filter(selection).Define("m_corr",tofill) \
                                             .Define("central_w_notnan", "TMath::IsNaN(central_w) ? 0.0 : central_w") \
                                             .Histo1D(model[0], "m_corr","central_w_notnan") 
          histos[var].Scale(1.0 / central_av) 
 

        elif sf_weights:
          # this only happens for data
          # sf_weight is a smart string
          print("YES I APPLY SF WEIGHTS!!!!!!")
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
            
          histos[var] = rdf.Filter(selection) \
                           .Define("central_w_notnan", "TMath::IsNaN(central_w) ? 0.0 : central_w") \
                           .Histo1D(model[0], tofill,"central_w_notnan") 
          histos[var].Scale(1.0 / central_av) 
 
        elif (sf_weights):
          # this only happens for data
          # sf_weight is a smart string
          print("YES I APPLY SF WEIGHTS!!!!!!")
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


def createBinnedPlots(splitter, regions, controlPlotsHighMass = None, controlPlotsLeftSideband = None, controlPlotsRightSideband = None, controlPlotsComplete = None, controlPlotsCustom = None  ):
  

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
  kk_wrong_incl     = f" && (k1_charge*k2_charge > 0)" #inclusive

  ################################################
  # Extract sigma from fit to MC DsMu signal to  #
  # define signal region.                        #
  ################################################

  # define signal region and sidebands :)
  mlow,   mhigh,   mlow2,   mhigh2,   mlow3,   mhigh3,   signalRegion,   anti_signalRegion,   leftSB,   rightSB   = getRegions(sigma)

  ################################
  # Fit the Ds mass for signflip #
  ################################

  # We fit the Ds mass on the most general region (just apply the baseline selection) 
  # and get form it the ratio of yields.
  # Later, if we continue with the binning loop, we can just apply this ratio, the selection
  # is automatically applied since the to-be-plotted histograms already care the selection  in its belly

  selection_massfit = baseline + right_sign 
  selec_massfit     = selections(selection_massfit)

  ## Signal and Hb (Hammer variations are never needed for the massfit)
  selec_M_DsMu            = createHistos(selec_massfit.dsMu       ,    rdfSig       , gen = False, hammer_central = hammer_central, hammer_sys = False ,sig = "dsmu"     )
  selec_M_DsMu_woHammer   = createHistos(selec_massfit.dsMu       ,    rdfSig       , gen = False)
  selec_M_DsTau           = createHistos(selec_massfit.dsTau      ,    rdfSig       , gen = False, hammer_central = hammer_central, hammer_sys = False ,sig = "dstau"    )
  selec_M_DsStarMu        = createHistos(selec_massfit.dsStarMu   ,    rdfSig       , gen = False, hammer_central = hammer_central, hammer_sys = False ,sig = "dsstarmu" )
  selec_M_DsStarTau       = createHistos(selec_massfit.dsStarTau  ,    rdfSig       , gen = False, hammer_central = hammer_central, hammer_sys = False ,sig = "dsstartau")
  selec_M_Hb              = createHistos(selec_massfit.hb         ,    rdfHb        , gen = False)
  selec_M_Mu_in_Hb        = createHistos(selec_massfit.dsMu       ,    rdfHb        , gen = False)
  selec_M_Data            = createHistos(selec_massfit.bare       ,    rdfData      , gen = False)

  selec_M_DsTau_blind     = { key: selec_M_DsTau[key].Clone()      for key in selec_M_DsTau.keys()     }
  selec_M_DsStarTau_blind = { key: selec_M_DsStarTau[key].Clone()  for key in selec_M_DsStarTau.keys() }

  #selec_M_DsTau["phiPi_m"].Scale(0.3)
  #selec_M_DsStarTau["phiPi_m"].Scale(0.3)

  #since the dsmu signal in Hb is not hammered, we calc. the scale for unhammered dsmu and hb. The effects should cancel! 
  hb_ratio_massfit = selec_M_Hb["phiPi_m"].Integral() / selec_M_Mu_in_Hb["phiPi_m"].Integral()
  hb_scale_massfit = selec_M_DsMu_woHammer["phiPi_m"].Integral() * hb_ratio_massfit 
   
 
  for key in selec_M_DsTau_blind.keys()    : selec_M_DsTau_blind[key]    .Scale(blind_scalar)    
  for key in selec_M_DsStarTau_blind.keys(): selec_M_DsStarTau_blind[key].Scale(blind_vector)    


  selec_M_Data_sf_pimu = createHistos(  baseline + data_selec + pimu_wrong    ,    rdfData       , gen = False, sf_weights = sf_weights , data = True, massfit = True)
  selec_M_Data_sf_kk   = createHistos(  baseline + data_selec + kk_wrong_incl ,    rdfData       , gen = False, sf_weights = sf_weights , data = True, massfit = True)

  # get the signflip scale by fitting the ds mass peak of sf data against hb + signal (called hRest)
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
  global scale_bkg,        scale_kk,       scale_pimu,       scale_n
  global scale_bkg_blind,  scale_kk_blind, scale_pimu_blind, scale_n_blind
 
  kk_prefit   = selec_M_Data_sf_kk  ["phiPi_m"].Clone().Integral()
  pimu_prefit = selec_M_Data_sf_pimu["phiPi_m"].Clone().Integral() 
  rest_prefit = hRest                          .Clone().Integral() 

  print("=============> Fit mass again to get signflip ratios")
  
  #if args.findcut: 

  #  bkg_postfit, rest_postfit       , abc = getSignflipRatioTest(selec_M_Data_sf_both["phiPi_m"].Clone(),  hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_m", 1.91, 2.028, findcut = False, cut = args.cut)


  hRest       = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                  [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                   selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                   selec_M_DsTau           ["phiPi_m"].Clone()  , 
                                   selec_M_DsStarTau       ["phiPi_m"].Clone()] ) 
 


  kk_postfit, pimu_postfit, rest_postfit       , abc = getSignflipRatio(selec_M_Data_sf_kk["phiPi_m"].Clone(),  selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_m", 1.91, 2.028)


  scale_kk         = kk_postfit         / kk_prefit 
  #scale_kk_blind   = kk_blind_postfit   / kk_prefit 
  scale_pimu       = pimu_postfit       / pimu_prefit 
  #scale_pimu_blind = pimu_blind_postfit / pimu_prefit 
  scale_rest       = rest_postfit       / rest_prefit 
  #scale_rest_blind = rest_blind_postfit / rest_prefit 


  # adapt selections for special regions

  if controlPlotsComplete:
    #just extend to the complete mass range to see the sidebands
    signalRegion = ""

  if controlPlotsCustom:
    #just extend with a custom selection 
    signalRegion += " && (disc_is_negative == 0)"


  if controlPlotsHighMass:
    # if we want to do control plots, we want to plot in the high mass region 
    # the splitter must be dsMu mass and the regions is [bs mass, 8]      
    splitter     = "dsMu_m"
    regions      = [[bsMass_,8]] #[[0,bsMass_],[bsMass_, 8]]
    low_mass     = "" # remove since we are now in the high mass region
    signalRegion = "" # remove to keep as much background as possible 
 
  if controlPlotsRightSideband:
    # if we want to do control plots, we want to plot in the sideband region 
    # the splitter must be phiPi mass and the regions leftSB and rightSB      
    splitter     = "phiPi_m"
    regions      = [[mhigh2,mhigh3]]
    print("Doing control plots, change regions to:", regions)
    signalRegion = "" # remove since we are now in the sidebands 
    low_mass     = "" # remove to keep as much background as possible  

  if controlPlotsLeftSideband:
    # if we want to do control plots, we want to plot in the sideband region 
    # the splitter must be phiPi mass and the regions leftSB and rightSB      
    splitter     = "phiPi_m"
    regions      = [[mlow3,mlow2]]
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
      #else it is a string, (used for 3D), f.e. region = "(q2_coll < 6) && (class == 0)"
      baseline_region   = baseline + f" && {region}" 
    

    # Define selections which hold all important selections next to baseline  
    selection        = baseline_region + right_sign + low_mass #for all other variables
    selection        += score_cut  
    data_selection   = baseline_region + data_selec + score_cut + pimu_wrong      + low_mass 

    ########################################
    # write normalization info into file   #
    ########################################

    norm_params = {
    "scale_pimu"      : scale_pimu,
    "scale_rest"      : scale_rest,
    "hb_ratio_massfit": hb_ratio_massfit
    }

    with open(f"{toSave_plots}/normalization_parameters.json", "w") as f:
      json.dump(norm_params,f)


    ###############################
    # Histograms in signal region #
    # NEW: submit to the queue    #
    ###############################

    #write a .sh file which can be submitted to the queue

    #top of the sh file
    core  = f"#!/bin/bash \n"
    core += f"eval \"$(conda shell.bash hook)\" \n"
    core += f"conda activate /work/pahwagne/environments/tf \n"

    #python command line
    python_command =  f"python plot_single_channel_and_bin.py" 
    python_command += f" --nn=\"{args.nn}\" --hammer=\"{args.hammer}\" --hammer_sys=\"{args.hammer_sys}\" --trigger=\"{args.trigger}\" --bdt=\"{args.bdt}\" "
    python_command += f" --bin={i} --region=\"{region}\" --toSave_plots=\"{toSave_plots}\" --selection=\"{selection}\" --data_selection=\"{data_selection}\" --split={split} "
    #python_command += " --scale_pimu={scale_pimu} --scale_rest={scale_rest} --hb_ratio_massfit={hb_ratio_massfit} "

    #optinal arguments
    if args.debug:   python_command += f" --debug"
    if args.control: python_command += f" --control={args.control}"
    if args.cut:     python_command += f" --cut={args.cut}" 

    #submitter command line
    sh_command = f"sbatch -p short "

    with open( toSave_plots + f"/submitter_DsMu_ch{i}.sh"         , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsMu ")
    os.system("chmod +x " + toSave_plots + f"/submitter_DsMu_ch{i}.sh" )
    #os.system(sh_command + toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  
    print    ( toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  
    os.system( toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsMu_woHammer_ch{i}.sh", "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsMu_woHammer")
    #os.system(sh_command + toSave_plots + f"/submitter_DsMu_woHammer_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsTau_ch{i}.sh"        , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsTau  ")
    #os.system(sh_command + toSave_plots + f"/submitter_DsTau_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsStarMu_ch{i}.sh"     , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsStarMu  ")
    #os.system(sh_command + toSave_plots + f"/submitter_DsStarMu_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsStarTau_ch{i}.sh"    , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsStarTau  ")
    #os.system(sh_command + toSave_plots + f"/submitter_DsStarTau_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Hb_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb")
    #os.system(sh_command + toSave_plots + f"/submitter_Hb_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Mu_in_Hb_ch{i}.sh"     , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Mu_in_Hb")
    #os.system(sh_command + toSave_plots + f"/submitter_Mu_in_Hb_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Data_ch{i}.sh"         , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Data")
    #os.system(sh_command + toSave_plots + f"/submitter_Data_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_sf_kk_ch{i}.sh"        , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Data_sf_pimu")
    #os.system(sh_command + toSave_plots + f"/submitter_sf_kk_ch{i}.sh"    )  
   
if control == "highmass":
  createBinnedPlots(split,binning, controlPlotsHighMass = True) 
elif control == "leftsb":
  createBinnedPlots(split,binning, controlPlotsLeftSideband = True) 
elif control == "rightsb":
  createBinnedPlots(split,binning, controlPlotsRightSideband = True) 
elif control == "complete":
  createBinnedPlots(split,binning, controlPlotsComplete = True) 
elif control == "custom":
  createBinnedPlots(split,binning, controlPlotsCustom = True)
else:
  createBinnedPlots(split,binning) 


