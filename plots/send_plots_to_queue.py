import glob
import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
import json

ROOT.ROOT.EnableImplicitMT(8) 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio, getSignflipRatioTest, fitAnotherVar
from helper import * 
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels, special_models, special_models_q2_coll,special_models_e_star_lhcb_alt
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
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
parser.add_argument("--bdt2",        required = True,     help = "Specify 'true' or 'false' to add bdt weights2") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
parser.add_argument("--control",                          help = "If given, run control plots, either 'highmass', 'leftsb', 'rightsb' or 'complete' or 'custom' ") 
parser.add_argument("--cut",                              help = "Cut on the discriminator score 5. ") 
parser.add_argument("--sel",                              help = "specify baseline selection (choose one in helper.py)") 
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

if args.bdt not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt , please use 'true' or 'false' (all lowercase!)")
else: bdt = (args.bdt == "true")

if args.bdt2 not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt2 , please use 'true' or 'false' (all lowercase!)")
else: bdt2 = (args.bdt2 == "true")

if args.sel not in baselines.keys() :
  raise ValueError ("Error: Not a valid key for --sel, please use one defined in helper.py")
else: sel = args.sel

if args.control and args.control not in ["highmass","leftsb","rightsb", "complete", "custom"]:
    raise ValueError ("Error: Not a valid key for --control, please use 'highmass', 'sb' or 'complete' or 'custom' ")
control = args.control


#optional arguments
#if args.cut  : score_cut = f" && (score5 <= {args.cut}) && (score1 > 0.1) && (q2_coll > 0) && (q2_coll < 12) && (score1 < 0.25)"
if args.cut and not pastNN : raise ValueError ("Error: Cannot interpret cut before NN")
elif args.cut and pastNN   : score_cut = f" && (score5 <= {args.cut} ) "
else                       : score_cut = ""

if args.debug: debug = 1000000
else         : debug = None

#update pastNN models
if pastNN: 
  models.update(pastNN_models)

#loading systematics
if hammer_sys:
  print(f"====> Adding the following systematics: {systematics_scalar} and {systematics_vector}") #defined in helper.py
  sys_scalar = systematics_scalar
  sys_vector = systematics_vector

with open("/work/pahwagne/RDsTools/hammercpp/development_branch/weights/20_10_2025_09_57_12/average_weights.yaml","r") as f:
  averages = yaml.safe_load(f)

# remove nan whammer weights:
hammer_str = ""

if hammer_central:
  hammer_str += "&& (TMath::Finite(central_w) > 0 )"

if hammer_sys:
  for sys in systematics_vector:
    hammer_str += f" && (TMath::Finite({sys}_up) > 0   )"
    hammer_str += f" && (TMath::Finite({sys}_down) > 0 )"


# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

# get date and time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S") 

#specify chan:variable for each chan to only plot this variable per chan
fitted_vars = {} #leave empty if you want all

#where to save plots and datacards
toSave_plots        = f"/work/pahwagne/RDsTools/plots/cmsplots_binned/{dt}/"
toSave_cards        = f"/work/pahwagne/RDsTools/fit/datacards_binned/{dt}"
shapes_folder       = f"/work/pahwagne/RDsTools/fit/shapes_binned/{dt}/"
shapes_folder_blind = f"/work/pahwagne/RDsTools/fit/shapes_binned/{dt}/blind/"

if not os.path.exists(toSave_plots): 
  os.makedirs(toSave_plots)
  os.makedirs(toSave_plots + "/log")  
  os.makedirs(toSave_plots + "/err")  

if not os.path.exists(shapes_folder): 
  os.makedirs(shapes_folder)

if not os.path.exists(shapes_folder_blind): 
  os.makedirs(shapes_folder_blind)

#also freeze sigma here!
sigma = 0.009
baseline = baselines[sel]

# Summary of settings
with open( toSave_plots + f"/info.txt", "a") as f:
  f.write(f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} on baseline {sel} on nn {nn_model}\n")
  print  (f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} on baseline {sel} on nn {nn_model}\n")


#############
# Split MC  #
#############

#if trigger == "mu7":
#  data_selec   = " && (mu7_ip4 == 1)"
#  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
#else:
#  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
#  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "

data_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#data_selec = " && ((mu7_ip4 == 1)||(mu9_ip6==1))"
#data_selec = " && ((mu7_ip4 == 1))"
mc_selec   = " "

# set the to be splitter variable and binning
split   = "cosMuW_lhcb_alt"
#split = "score2"

binning = [[-99,99]]
#binning = [[0,0.7],[0.7,1]]

################
## 3D binning ##
################

binning_str = []

sr_string = " && ((phiPi_m >= 1.94134) && (phiPi_m <= 1.99534 ) )"
sb_string = " && ((phiPi_m < 1.94 ) || (phiPi_m > 1.995  ) )"

#used later
mass_sb     = "    ((phiPi_m <  1.92) || (phiPi_m >  2.02))" 
mass_center = " && ((phiPi_m >= 1.92) && (phiPi_m <= 2.02))"

#comb control region
#binning_str.append( mass_sb                                                                                                          ); fitted_vars[0]  = "phiPi_m";
##
###dsmu control region
##binning_str.append("(score2 > 0.3)                                       && (q2_coll >= 0)  && (q2_coll < 2) "  + mass_center); fitted_vars[1]  = "score2";
##binning_str.append("(score2 > 0.3)                                       && (q2_coll >= 2)  && (q2_coll < 4) "  + mass_center); fitted_vars[2]  = "score2";
##binning_str.append("(score2 > 0.3)                                       && (q2_coll >= 4)  && (q2_coll < 6) "  + mass_center); fitted_vars[3]  = "score2";
##binning_str.append("(score2 > 0.3)                                       && (q2_coll >= 6)  && (q2_coll < 8) "  + mass_center); fitted_vars[4]  = "score2";
##binning_str.append("(score2 > 0.3)                                       && (q2_coll >= 8)  && (q2_coll < 10)"  + mass_center); fitted_vars[5]  = "score2";
#binning_str.append("(score2 > 0.6)                                       "  + mass_center); fitted_vars[1]  = "q2_coll";
##
##ds*mu center region
##binning_str.append("(score2 < 0.3) && (score3 > 0.4)                     && (q2_coll >= 0)  && (q2_coll < 2) "  + mass_center); fitted_vars[6]  = "score2";
##binning_str.append("(score2 < 0.3) && (score3 > 0.4)                     && (q2_coll >= 2)  && (q2_coll < 4) "  + mass_center); fitted_vars[7]  = "score2";
##binning_str.append("(score2 < 0.3) && (score3 > 0.4)                     && (q2_coll >= 4)  && (q2_coll < 6) "  + mass_center); fitted_vars[8]  = "score2";
##binning_str.append("(score2 < 0.3) && (score3 > 0.4)                     && (q2_coll >= 6)  && (q2_coll < 8)"  + mass_center); fitted_vars[9]   = "score2";
##binning_str.append("(score2 < 0.3) && (score3 > 0.4)                     && (q2_coll >= 8)  && (q2_coll < 10)"  + mass_center); fitted_vars[10] = "score2";
#binning_str.append("(score2 < 0.6) && (score3 > 0.4)                     "  + mass_center); fitted_vars[2]  = "q2_coll";
##
##hb control region
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 > 0.35)                                              "  + mass_center); fitted_vars[3] = "score1";
#
###sr
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 0)  && (q2_coll < 2) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[4] = "score1";
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 0)  && (q2_coll < 4) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[5] = "score1";
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 4)  && (q2_coll < 5) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[6] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 0)  && (q2_coll < 6) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[4] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 6)  && (q2_coll < 7) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[5] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 7)  && (q2_coll < 8) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[6] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 8)  && (q2_coll < 9) "  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[7] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 9)  && (q2_coll < 10)"  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[8] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35)  && (q2_coll >= 10) && (q2_coll < 12)"  + " && (phiPi_m > 1.96) && (phiPi_m < 1.976) "); fitted_vars[9] = "score1";
##
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 0)  && (q2_coll < 2) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[13] = "score1";
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 2)  && (q2_coll < 3) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[14] = "score1";
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 3)  && (q2_coll < 4) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[15] = "score1";
###binning_str.append("(score2 < 0.3) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 4)  && (q2_coll < 5) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[16] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 0)  && (q2_coll < 6) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[10] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 6)  && (q2_coll < 7) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[11] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 7)  && (q2_coll < 8) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[12] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 8)  && (q2_coll < 9) " + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[13] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 9)  && (q2_coll < 10)" + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[14] = "score1";
#binning_str.append("(score2 < 0.6) && (score3 < 0.4) && (score4 < 0.35) && (q2_coll >= 10) && (q2_coll < 12)" + " && (((phiPi_m > 1.92) && (phiPi_m < 1.96)) || ((phiPi_m > 1.976) && (phiPi_m < 2.02))) "); fitted_vars[15] = "score1";

##########################################

for b in binning:
  
  # b is a list
  #binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" + sr_string)  
  #binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" + sb_string)  
  binning_str.append( f"( {split} > {b[0]} ) && ( {split} < {b[1]} )" )  

##overwrite
binning = binning_str

if ((len(fitted_vars) != 0) and (len(fitted_vars) != len(binning))): print("Abort! Bin:var dict doesn't match binning!"); sys.exit()


##########################
# write info into file   #
##########################

with open(f"{toSave_plots}/fitted_vars.json", "w") as f:
  json.dump(fitted_vars,f)

##############################
# Load chain into RDataFrame #
##############################

def getRdf(files, path, debug = None, bph_part = None, isData = None):

  # prepare the chain
  chain = ROOT.TChain("tree")

  if debug: 

    files_str = f"{path}/{files[0]}/*1*.root"
    n = chain.Add(files_str)

  else:    

    files_str = [f"{path}/{f}/*.root" for f in files]
    n = chain.Add(files_str)

  rdf = ROOT.RDataFrame("tree", files_str)
  print(f"====> Adding {n} files for plotting")

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
print(f"====> Start creating RDataFrames")

#prepare the defaults
files_sig  = sig_flatNanos 
files_hb   = hb_flatNanos
files_data = data_flatNanos 

#default: we always run on pre-skimmed files
path_sig  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"
path_hb   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"
path_data = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"

if (not pastNN and hammer_central):

  files_sig  = [sig_hammer_flatNano]
  path_sig   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/" 

if (pastNN and not bdt and not bdt2):

  files_sig  = [sig_pastNN ]
  files_hb   = [hb_pastNN  ]
  files_data = [data_pastNN]

  path_sig   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees_test/" 
  path_hb    = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees_test/" 
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees_test/" 


if (pastNN and bdt and not bdt2):

  file_data  = [bdt_data] 
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/" 

if (pastNN and bdt and bdt2):

  file_data  = [bdt_data_afternn] 
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/" 


if (pastNN and not bdt and bdt2):

  print("Illegal argument combination bdt {bdt} and bdt2 {bdt2} !")
  sys.exit()


print(f"====> Running on the following files: \n sig: {path_sig} \n hb: {path_hb} \n data: {path_data} \n")


chainSigSB, rdfSigSB     = getRdf( files_sig , path_sig , debug = debug)
chainSig,   rdfSig       = getRdf( files_sig , path_sig , debug = debug)
chainHb,    rdfHb        = getRdf( files_hb  , path_hb  , debug = debug)
_,          rdfHbHammer  = getRdf( files_hb  , path_hb  , debug = debug)
chainData,  rdfData      = getRdf( files_data, path_data, debug = debug)

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
  hb.Scale(hb_scale)

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

def createHistos(selection, rdf, data = False , variables = None, ff_central = False, ff_sys = False, sig = None, sf_weights = None, sf_weights2 = None, region = None, massfit = False):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"
  "If <variables> are given, only these variables from <models> are created" #f.e. for the phiPi mass or DsMu with different selection

  print(" ====> ff_sys: "        , ff_sys   )
  print(" ====> selection is: "  , selection)

  histos = {}

  if pastNN:           models.update(pastNN_models)
  if massfit == False: models.update(modelsSR)

  # define the final weight string

  total_w_str = ""

  if (sig and not ff_central):
    total_w_str = "trigger_sf"
  
  if (sig and ff_central):

    central_av = averages[ "central_w_" + sig]
    total_w_str = f"trigger_sf * central_w / {central_av}"
 
  if (data and sf_weights and not sf_weights2):
    total_w_str = "sf_weights"

  if (data and sf_weights and sf_weights2):
    total_w_str = "sf_weights * sf_weights2"

  if data:
    #for safety, when we plot data we dont apply any weights!!!
    total_w_str = ""

  print(f"====> Total applied weight is: {total_w_str}" )

  for var,model in models.items(): 

      ################################
      # Special binnings             #
      ################################

      if var == "score1" and (split == "class") and region != None:
        
        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")
  
        model = special_models[var + f"_bin{region}" ]

      if var == "score1" and (split == "q2_coll" or split == "q2_lhcb_alt") and region != None:
        
        #adapt the binning
        print(f"======> Adapt binning for {var} and region {region}")
  
        model = special_models_q2_coll[var + f"_bin{region}" ]                                     

      if var == "e_star_lhcb_alt" and (split == "q2_coll" or split == "q2_lhcb_alt") and region != None:
        
        #adapt the binning
        print(f"======> Adapt binning for {var} and region {region}")
  
        model = special_models_e_star_lhcb_alt[var + f"_bin{region}" ]                                     


      ##############################################
      # Fill only variables, if explicitly given!  #
      ##############################################

      if variables:
        if var not in variables: 
          print(f"Skip plotting variable {var}");     
          continue 


      ##############################################
      # Start filling                              #
      ##############################################

      print("Filling variable ", var)

      # Shift mass peak of MC by half permill
      #bool      #true for MC                    #true for data
      mass_expr   = f"(run==1) * 0.9995 * phiPi_m + (run!=1) * phiPi_m"
 
      if total_w_str == "":
        histos[var] = rdf.Filter(selection).Define("m_corr",mass_expr).Histo1D(model[0], var ) 
      else:
        histos[var] = rdf.Filter(selection).Define("m_corr",mass_expr).Define("total_w", total_w_str).Histo1D(model[0], var , "total_w") 
 
      if ff_sys:

        #this is only true for signals
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
          histos[var + "_" + s + "Up"]    = rdf.Filter(selection).Define("m_corr", mass_expr).Histo1D(model[0], var, s + "_up"   )
          histos[var + "_" + s + "Up"]    .Scale(1.0 / var_up_av)

          histos[var + "_" + s + "Down"]  = rdf.Filter(selection).Define("m_corr", mass_expr).Histo1D(model[0], var, s + "_down" )
          histos[var + "_" + s + "Down" ] .Scale(1.0 / var_down_av)
            

      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetBinContent(histos[var].GetMaximumBin()))
      
      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(2)

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

  selection_massfit = baseline + right_sign + low_mass
  selec_massfit     = selections(selection_massfit)

  print("---> before asking for selec_M")

  ## Signal and Hb (Hammer variations are never needed for the massfit)
  print("---> filling DsMu")
  selec_M_DsMu            = createHistos(selec_massfit.dsMu      + score_cut + hammer_str,    rdfSig       ,variables = ["phiPi_m"] , ff_central = hammer_central, sig = "dsmu"     )
  print("---> filling DsMu")
  selec_M_DsMu_woHammer   = createHistos(selec_massfit.dsMu      + score_cut + hammer_str,    rdfSig       ,variables = ["phiPi_m"] ,                              sig = "dsmu"     )
  print("---> filling DsTau")
  selec_M_DsTau           = createHistos(selec_massfit.dsTau     + score_cut + hammer_str,    rdfSig       ,variables = ["phiPi_m"] , ff_central = hammer_central, sig = "dstau"    )
  print("---> filling DsStarMu")
  selec_M_DsStarMu        = createHistos(selec_massfit.dsStarMu  + score_cut + hammer_str,    rdfSig       ,variables = ["phiPi_m"] , ff_central = hammer_central, sig = "dsstarmu" )
  print("---> filling DsStarTau")
  selec_M_DsStarTau       = createHistos(selec_massfit.dsStarTau + score_cut + hammer_str,    rdfSig       ,variables = ["phiPi_m"] , ff_central = hammer_central, sig = "dsstartau")


  print("---> filling Hb")
  selec_M_Hb              = createHistos(selec_massfit.hb        + score_cut             ,    rdfHb        ,variables = ["phiPi_m"]                                             )
  print("---> filling Hb")
  selec_M_Mu_in_Hb        = createHistos(selec_massfit.dsMu      + score_cut             ,    rdfHb        ,variables = ["phiPi_m"]                                             )
  print("---> filling Data")
  selec_M_Data            = createHistos(selec_massfit.bare      + score_cut             ,    rdfData      ,variables = ["phiPi_m"] , data = True                               )

  selec_M_DsTau_blind     = { key: selec_M_DsTau[key].Clone()      for key in selec_M_DsTau.keys()     }
  selec_M_DsStarTau_blind = { key: selec_M_DsStarTau[key].Clone()  for key in selec_M_DsStarTau.keys() }


  #selec_M_DsTau["phiPi_m"].Scale(0.3)
  #selec_M_DsStarTau["phiPi_m"].Scale(0.3)

  #since the dsmu signal in Hb is not hammered, we calc. the scale for unhammered dsmu and hb. The effects should cancel! 
  hb_ratio_massfit = selec_M_Hb["phiPi_m"].Integral() / selec_M_Mu_in_Hb["phiPi_m"].Integral()
  hb_scale_massfit = selec_M_DsMu["phiPi_m"].Integral() * hb_ratio_massfit / selec_M_Hb["phiPi_m"].Integral() 
   
 
  for key in selec_M_DsTau_blind.keys()    : selec_M_DsTau_blind[key]    .Scale(blind_scalar)    
  for key in selec_M_DsStarTau_blind.keys(): selec_M_DsStarTau_blind[key].Scale(blind_vector)    


  #no BDT2 correction here since we do this plot before the nn cut :D
  print("---> filling pimu flip")
  selec_M_Data_sf_pimu = createHistos(  baseline + data_selec + pimu_wrong    + score_cut + low_mass,    rdfData       , variables= ["phiPi_m"], sf_weights = bdt, sf_weights2 = bdt2 ) 

  pimu_prefit = selec_M_Data_sf_pimu["phiPi_m"].Clone().Integral() 


  # get the signflip scale by fitting the ds mass peak of sf data against hb + signal (called hRest)

  hRest_blind = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                 [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                  selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                  selec_M_DsTau_blind     ["phiPi_m"].Clone()  , 
                                  selec_M_DsStarTau_blind ["phiPi_m"].Clone()] ) 

  rest_prefit_blind = hRest_blind                    .Clone().Integral() 

  pimu_postfit_blind, rest_postfit_blind, abc = getSignflipRatio(selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest_blind      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_m", 1.91, 2.028)


  hRest       = prepareSignFlip(  selec_M_Hb              ["phiPi_m"].Clone()  , hb_scale_massfit,
                                  [selec_M_DsMu            ["phiPi_m"].Clone()  , 
                                   selec_M_DsStarMu        ["phiPi_m"].Clone()  , 
                                   selec_M_DsTau           ["phiPi_m"].Clone()  , 
                                   selec_M_DsStarTau       ["phiPi_m"].Clone()] ) 

  rest_prefit = hRest                          .Clone().Integral() 

  pimu_postfit, rest_postfit, abc = getSignflipRatio(selec_M_Data_sf_pimu["phiPi_m"].Clone(), hRest      ,selec_M_Data["phiPi_m"].Clone(), mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, "phiPi_m", 1.91, 2.028)

  # now get the yiel ratios (postfit/prefit) 
  global scale_bkg,        scale_pimu,       scale_n
  global scale_bkg_blind,  scale_pimu_blind, scale_n_blind


  scale_pimu       = pimu_postfit       / pimu_prefit 
  scale_pimu_blind = pimu_postfit_blind / pimu_prefit 

  scale_rest       = rest_postfit       / rest_prefit 
  scale_rest_blind = rest_postfit_blind / rest_prefit_blind 


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


    with open( toSave_plots + f"/info.txt", "a") as f:
      f.write( f"====> ch{i} belongs to region: {region} \n")
  
    

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
    "scale_pimu"            : scale_pimu,
    "scale_pimu_blind"      : scale_pimu_blind,
    "scale_rest"            : scale_rest,
    "scale_rest_blind"      : scale_rest_blind,

    "hb_ratio_massfit": hb_ratio_massfit,
    "hb_scale_massfit": hb_scale_massfit,
    #"trigger"         : trigger,
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
    #python_command += f" --nn=\"{args.nn}\" --hammer=\"{args.hammer}\" --hammer_sys=\"{args.hammer_sys}\" --trigger=\"{args.trigger}\" --bdt=\"{args.bdt}\" --bdt2=\"{args.bdt2}\" "
    python_command += f" --nn=\"{args.nn}\" --hammer=\"{args.hammer}\" --hammer_sys=\"{args.hammer_sys}\" --bdt=\"{args.bdt}\" --bdt2=\"{args.bdt2}\" "
    python_command += f" --bin={i} --region=\"{region}\" --toSave_plots=\"{toSave_plots}\" --selection=\"{selection}\" --data_selection=\"{data_selection}\" --split={split} "
    #python_command += " --scale_pimu={scale_pimu} --scale_rest={scale_rest} --hb_ratio_massfit={hb_ratio_massfit} "

    #optinal arguments
    if args.debug:   python_command += f" --debug"
    if args.control: python_command += f" --control={args.control}"
    if args.cut:     python_command += f" --cut={args.cut}" 

    #submitter command line
    #sh_command = f"sbatch -p short -o {toSave_plots}/log/log{i}.log -e {toSave_plots}/err/err{i}.err"
    sh_command = f"sbatch -p short "

    with open( toSave_plots + f"/submitter_DsMu_ch{i}.sh"         , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsMu ")
    os.system("chmod +x " + toSave_plots + f"/submitter_DsMu_ch{i}.sh" )
    os.system(sh_command + f"-o {toSave_plots}/log/log_DsMu{i}.log -e {toSave_plots}/err/err_DsMu{i}.err " + toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  
    #print(sh_command + f"-o {toSave_plots}/log/log_DsMu{i}.log -e {toSave_plots}/err/err_DsMu{i}.err " + toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  
    #print    ( toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  
    #os.system( toSave_plots + f"/submitter_DsMu_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsMu_woHammer_ch{i}.sh", "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsMu_woHammer")
    os.system(sh_command + f" -o {toSave_plots}/log/log_DsMu{i}.log -e {toSave_plots}/err/err_DsMu{i}.err " + toSave_plots + f"/submitter_DsMu_woHammer_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsTau_ch{i}.sh"        , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsTau  ")
    os.system(sh_command + f"-o {toSave_plots}/log/log_DsTau{i}.log -e {toSave_plots}/err/err_DsTau{i}.err " + toSave_plots + f"/submitter_DsTau_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsStarMu_ch{i}.sh"     , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsStarMu  ")
    os.system(sh_command + f"-o {toSave_plots}/log/log_DsStarMu{i}.log -e {toSave_plots}/err/err_DsStarMu{i}.err " + toSave_plots + f"/submitter_DsStarMu_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_DsStarTau_ch{i}.sh"    , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=DsStarTau  ")
    os.system(sh_command + f"-o {toSave_plots}/log/log_DsStarTau{i}.log -e {toSave_plots}/err/err_DsStarTau{i}.err " + toSave_plots + f"/submitter_DsStarTau_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Hb_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb{i}.log -e {toSave_plots}/err/err_Hb{i}.err " + toSave_plots + f"/submitter_Hb_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_fd_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_fd")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_fd{i}.log -e {toSave_plots}/err/err_Hb_fd{i}.err " + toSave_plots + f"/submitter_Hb_fd_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_dc_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_dc")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_dc{i}.log -e {toSave_plots}/err/err_Hb_dc{i}.err " + toSave_plots + f"/submitter_Hb_dc_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_others_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_others")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_others{i}.log -e {toSave_plots}/err/err_Hb_others{i}.err " + toSave_plots + f"/submitter_Hb_others_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_bs_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_bs")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_bs{i}.log -e {toSave_plots}/err/err_Hb_bs{i}.err " + toSave_plots + f"/submitter_Hb_bs_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_bpm_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_bpm")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_bpm{i}.log -e {toSave_plots}/err/err_Hb_bpm{i}.err " + toSave_plots + f"/submitter_Hb_bpm_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_b0_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_b0")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_b0{i}.log -e {toSave_plots}/err/err_Hb_b0{i}.err " + toSave_plots + f"/submitter_Hb_b0_ch{i}.sh"    )  

    with open( toSave_plots + f"/submitter_Hb_lambdab_ch{i}.sh"           , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Hb_lambdab")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Hb_lambdab{i}.log -e {toSave_plots}/err/err_Hb_lambdab{i}.err " + toSave_plots + f"/submitter_Hb_lambdab_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Mu_in_Hb_ch{i}.sh"     , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Mu_in_Hb")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Mu_in_Hb{i}.log -e {toSave_plots}/err/err_Mu_in_Hb{i}.err " + toSave_plots + f"/submitter_Mu_in_Hb_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_Data_ch{i}.sh"         , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Data")
    os.system(sh_command + f"-o {toSave_plots}/log/log_Data{i}.log -e {toSave_plots}/err/err_Data{i}.err " + toSave_plots + f"/submitter_Data_ch{i}.sh"    )  


    with open( toSave_plots + f"/submitter_sf_pimu_ch{i}.sh"        , "w") as f:
      f.write(core)
      f.write(f" {python_command} --channel=Data_sf_pimu")
    os.system(sh_command + f"-o {toSave_plots}/log/log_pimu{i}.log -e {toSave_plots}/err/err_pimu{i}.err " + toSave_plots + f"/submitter_sf_pimu_ch{i}.sh"    )  
   
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


