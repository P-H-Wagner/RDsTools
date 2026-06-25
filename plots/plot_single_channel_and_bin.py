import array
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
from signflip  import getSignflipRatio, getSignflipRatioTest, fitAnotherVar
from helper import * 
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels, special_models_score2, special_models_bdt
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
parser.add_argument("--nn",          required = True,     help = "Specify 'before' or 'past' neural network plots ") 
parser.add_argument("--hammer",      required = True,     help = "Specify 'true' or 'false' to apply hammer weights") 
parser.add_argument("--hammer_sys",  required = True,     help = "Specify 'true' or 'false' to save weight variation shapes") 
parser.add_argument("--bs_tau",      required = True,     help = "Specify 'true' or 'false' to apply bs lifetime weights") 
parser.add_argument("--bs_tau_sys",  required = True,     help = "Specify 'true' or 'false' to save bs lifetime weight variation shapes") 
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
parser.add_argument("--bdt2",         required = True,     help = "Specify 'true' or 'false' to add bdt2 weights") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
#parser.add_argument("--findcut",     action='store_true', help = "If given, we run thecut scan") 

parser.add_argument("--bin",   type=int, required = True,     help = "") 
parser.add_argument("--region",          required = True,     help = "") 
parser.add_argument("--toSave_plots",    required = True,     help = "") 
parser.add_argument("--selection",       required = True,     help = "") 
parser.add_argument("--data_selection",  required = True,     help = "") 
parser.add_argument("--channel",         required = True,     help = "") 
parser.add_argument("--split",           required = True,     help = "") 
parser.add_argument("--parts",           required = True,     help = "specify which parts to plot, comma separated, f.e.: 1,2,3")
parser.add_argument("--model",                                help = "overwrite model in helper file")



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

if args.bs_tau not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bs_tau, please use 'true' or 'false' (all lowercase!)")
else: bs_tau_central = (args.bs_tau == "true")

if args.bs_tau_sys not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --sys, please use 'true' or 'false' (all lowercase!)")
else: bs_tau_sys = (args.bs_tau_sys == "true")


if args.bdt not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt , please use 'true' or 'false' (all lowercase!)")
else: bdt = (args.bdt == "true")

if args.bdt2 not in ["true", "false"]:
  raise ValueError ("Error: Not a valid key for --bdt2 , please use 'true' or 'false' (all lowercase!)")
else: bdt2 = (args.bdt2 == "true")


if args.parts :
  parts = args.parts.split(",")
  parts = [int(x) for x in parts]
else:
  parts = None

if args.model:
  print(f"====> Overwriting nn modek, using: {args.model} instead")
  nn_model = args.model


if args.debug: debug = 50000
else         : debug = None




if pastNN: models.update(pastNN_models)

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
    hammer_str += f" && (TMath::Finite({sys}_up) > 0 )"
    hammer_str += f" && (TMath::Finite({sys}_down) > 0 )"


# disable title and stats and displaying
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!

# get date and time
now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S") 

#where to save plots and datacards

# Summary of settings
with open( args.toSave_plots + f"/info.txt", "a") as f:
  f.write(f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} with bs lt: {bs_tau_central} and bs lt sys: {bs_tau_sys} \n")
  print  (f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} with bs lt: {bs_tau_central} and bs lt sys: {bs_tau_sys} \n")
  

#############
# Split MC  #
#############

#if trigger == "mu7":
#  data_selec   = " && (mu7_ip4 == 1)"
#  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
#else:
#  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
#  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "

#data_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8p5_3p5)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu10p5_3p5)||(mu12_ip6))"
#data_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#data_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#data_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#data_selec = " && ((mu7_ip4 ==0) && (mu9_ip6==1))"
#data_selec = " && ((mu7_ip4 != 1) && (mu9_ip6!=1) && (mu12_ip6 == 1))"
#data_selec = " && ((mu7_ip4 == 1))"
#data_selec = " && ((mu8_ip3 == 1))"
#data_selec = " && ((mu8_ip5 == 1))"
#data_selec = " && ((mu8_ip6 == 1))"
#data_selec = " && ((mu9_ip4 == 1))"
#data_selec = " && ((mu9_ip5 == 1))"
#data_selec = " && ((mu9_ip6 == 1))"
#data_selec = " && ((mu12_ip6 == 1))"
data_selec = " "


#mc_selec  = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8p5_3p5)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu10p5_3p5)||(mu12_ip6))"
#mc_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "
#mc_selec = " && ((mu7_ip4 == 1))"
#mc_selec = " && ((mu8_ip3 == 1))"
#mc_selec = " && ((mu8_ip5 == 1))"
#mc_selec = " && ((mu8_ip6 == 1))"
#mc_selec = " && ((mu9_ip4 == 1))"
#mc_selec = " && ((mu9_ip5 == 1))"
#mc_selec = " && ((mu9_ip6 == 1))"
#mc_selec  = " && ((mu12_ip6 == 1))"
mc_selec  = " "

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
    for f in files_str: n = chain.Add(f)


  rdf = ROOT.RDataFrame("tree", files_str)
  print(f"====> Adding {n} files for plotting:")
  print(files_str)

  return (chain,rdf)


########################################
# Assign selections depending on the   #
# MC signal ID                         #
########################################

class selections:

  def __init__(self, selec):


    self.hb =        selec + mc_selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1)  && (gen_same_mother == 1) "
    self.hb_dc =     selec + mc_selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1)  && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 != 7) && (gen_sig < 500)" 
    self.hb_fd =     selec + mc_selec + " && (gen_sig != 0) && (gen_sig != 1) && (gen_sig != 10) && (gen_sig != 11) & (gen_match_success ==1)  && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 == 7)"
    self.hb_others = selec + mc_selec + " && (gen_match_success ==1)  && (gen_same_mother == 1) && (gen_sig >= 500)"

    self.hb_bs    =  selec + mc_selec + " && (300 <= gen_sig) && (gen_sig < 400) && (gen_same_mother == 1) "
    self.hb_bs_fd =  selec + mc_selec + " && (300 <= gen_sig) && (gen_sig < 400) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 == 7)"
    self.hb_bs_dc =  selec + mc_selec + " && (300 <= gen_sig) && (gen_sig < 400) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 != 7)"

    self.hb_b0    =     selec + mc_selec + " && (200 <= gen_sig) && (gen_sig < 300) && (gen_same_mother == 1) "
    self.hb_b0_fd =     selec + mc_selec + " && (200 <= gen_sig) && (gen_sig < 300) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 == 7)"
    self.hb_b0_dc =     selec + mc_selec + " && (200 <= gen_sig) && (gen_sig < 300) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 != 7)"


    self.hb_bpm    =    selec + mc_selec + " && (100 <= gen_sig) && (gen_sig < 200) && (gen_same_mother == 1) "
    self.hb_bpm_fd =    selec + mc_selec + " && (100 <= gen_sig) && (gen_sig < 200) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 == 7)"
    self.hb_bpm_dc =    selec + mc_selec + " && (100 <= gen_sig) && (gen_sig < 200) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 != 7)"

    self.hb_lambdab    = selec + mc_selec + " && (400 <= gen_sig) && (gen_sig < 500) && (gen_same_mother == 1) "
    self.hb_lambdab_fd = selec + mc_selec + " && (400 <= gen_sig) && (gen_sig < 500) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 == 7)"
    self.hb_lambdab_dc = selec + mc_selec + " && (400 <= gen_sig) && (gen_sig < 500) && (gen_same_mother == 1) && (static_cast<int>(gen_sig) % 10 != 7)"

    self.dsMu =      selec + mc_selec + " && (gen_sig == 0)"   
    self.dsTau =     selec + mc_selec + " && (gen_sig == 1)"   
    self.dsStarMu =  selec + mc_selec + " && (gen_sig == 10)"  
    self.dsStarTau = selec + mc_selec + " && (gen_sig == 11)"  

    self.bare =      selec + data_selec


#########################################
## INPUT SAMPLES                       ##
#########################################

# Create rdf from tree
print(f" ===> Start creating RDataFrames")


#prepare the defaults
files_sig  = sig_flatNanos
files_hb   = hb_flatNanos
files_data = data_flatNanos

#default: we always run on pre-skimmed files
path_sig  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"
path_hb   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"
path_data = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/"

#update sig
if hammer_central:

  files_sig  = [sig_hammer_flatNano]
  path_sig   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/"

#update sig and data (includes hammer weights and bdt weights anyways)
if pastNN:

  files_sig  = [f"sig_{nn_model}"]
  files_hb   = [f"hb_{nn_model}" ]
  if parts:
    files_data = [f"data_bph{part}_{nn_model}" for part in parts]
  else:
    files_data = [f"data_{nn_model}" ]

  path_sig   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/"
  path_hb    = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/"
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/"

#update data
if (bdt and not bdt2):

  if parts:  
    files_data = [f"{bdt_model}/data_bph{part}_{bdt_model}/" for part in parts]
  else:
    files_data = [f"{bdt_model}" ]

  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/" 


#update data again (bdt2 only exists with nn!)
if (pastNN and bdt and bdt2):

  files_data = [bdt_data_afternn]
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/"


if (pastNN and not bdt and bdt2) or (not pastNN and bdt2):

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


def createHistos(selection, rdf, data = False , variables = None, ff_central = False, ff_sys = False, lt_central = False, lt_sys = False, mc = None, sf_weights = None, sf_weights2 = None, region = None, massfit = False, comb = None, combSys = False):

  "Creates histograms of all histograms in <models> (prepared in histModels.py) with the <selection>"
  "If <variables> are given, only these variables from <models> are created" #f.e. for the phiPi mass or DsMu with different selection

  print(" ====> ff_sys: "        , ff_sys   )
  print(" ====> selection is: "  , selection)

  histos = {}

  if pastNN:           models.update(pastNN_models)
  if massfit == False: models.update(modelsSR)

  # define the final weight string

  total_w_str = ""

  if (mc and not ff_central):
    total_w_str = "trigger_sf"
    #total_w_str = ""

  if (mc and ff_central):

    central_av = averages[ "central_w_" + mc]
    total_w_str = f"trigger_sf * central_w / {central_av}"
    #total_w_str = f"central_w / {central_av}"

  if (sf_weights and not sf_weights2):
    total_w_str = "sf_weights"

  if (sf_weights and sf_weights2):
    total_w_str = "sf_weights * sf_weights2"

  #apply bs lifetime weights for signal (mc given, but not hb)
  if (mc and mc != "hb" and lt_central):

     if total_w_str != "": total_w_str += " * w_bs_tau"
     else                : total_w_str += "   w_bs_tau"


  print(f"====> Total applied weight is: {total_w_str}" )

  ################################################################
  # Filter rdf and modify rdf once and for all for all variables #
  ################################################################

  # Shift mass peak of MC by half permill --> weights that change along x

  # This is not a per-event weight! Just a change along x, only applied for the mass itself
  #bool          #true for MC                            #true for data
  mass_shifted         = f"(run==1) * 0.9995 * phiPi_m         + (run!=1) * phiPi_m"
  mass_shifted_smeared = f"(run==1) * 0.9995 * smear * phiPi_m + (run!=1) * phiPi_m"
  #mass_expr    = f"(run==1) * 0.9995 * phiPi_m + (run!=1) * phiPi_m"

  # Tilt mass shape of combinatorial --> weights that change along y
  # This is a per event weights! Change along y, needs to be applied for all variables
  if comb:

    y0 = n0_exp * np.exp(tau * m0)
    y1 = n0_exp * np.exp(tau * m1)

    y0_new = y0 - tilt * y0
    y1_new = y1 + tilt * y1

    #define a linear function, tilted
    a  = (y1_new - y0_new) / (m1 - m0)
    b  = y0_new - a * m0

    if total_w_str != "":
      total_w_str  += f" * (( {a} * phiPi_m + {b} ) / ( {n0_exp} * std::exp({tau} * phiPi_m )))"
    else:
      total_w_str  += f"   (( {a} * phiPi_m + {b} ) / ( {n0_exp} * std::exp({tau} * phiPi_m )))"



    print(f"====> tilting combinatorial by applying weight {total_w_str}")

  #correct dxy
  dxy_err_expr = f"(run==1) * 1.0 * dxy_mu_err_pv + (run!=1) * dxy_mu_err_pv"
  dxy_sig_expr = f" dxy_mu_pv / dxy_mu_err_pv_corr "

  #re-order k1 and k2 after pt
  k1_expr      = f" k1_pt * (k1_pt > k2_pt) + k2_pt * (k1_pt <= k2_pt) "
  k2_expr      = f" k2_pt * (k1_pt > k2_pt) + k1_pt * (k1_pt <= k2_pt) "

  rdf_filt = rdf\
             .Define("m_shifted",mass_shifted)\
             .Filter(selection)\
             .Define("smear","gRandom->Gaus(1.0, 0.001)")\
             .Define("m_corr"   ,mass_shifted_smeared)\
             .Define("dxy_mu_err_pv_corr",dxy_err_expr)\
             .Define("dxy_mu_sig_pv_corr",dxy_sig_expr)\
             .Define("k1_corr",k1_expr)\
             .Define("k2_corr",k2_expr)\



  for var,model in models.items():


      #in the case we have a special binning, only plot the corresponding variable for the region!
      if ((len(fitted_vars) != 0) and (var != fitted_vars[str(region)])): continue

      ################################
      # Special binnings             #
      ################################

      #if ("score" in var or var == "phiPi_m" or "q2" in var) and (args.split == "score2" ) and region != None:
      #  
      #  #adapt the binning
      #  print(f"======> Adapt binning for {var} and region {region}")
  
      #  model = special_models_score2[var + f"_bin{region}" ]  
      #  print(model)


      model = special_models_bdt[var + f"_bin{region}" ]
      

      ##############################################
      # Fill only variables, if explicitly given!  #
      ##############################################

      if variables:
        if var not in variables:
          print(f"Skip plotting variable {var}");
          continue

      ##############################################
      # Assign combinatorial systematic            #
      #############################################

      # this has to happen in the variable loop bc the histograms are per variable!

      if combSys:
  
        with open(f"/work/pahwagne/RDsTools/corrections/combinatorial/25_06_2026_13_41_29/combSys/comb_sys_weights_up_{var}_ch{region}.json"  ,"r") as f:
          comb_up = json.load(f) 
        with open(f"/work/pahwagne/RDsTools/corrections/combinatorial/25_06_2026_13_41_29/combSys/comb_sys_weights_down_{var}_ch{region}.json","r") as f:
          comb_down = json.load(f) 
        with open(f"/work/pahwagne/RDsTools/corrections/combinatorial/25_06_2026_13_41_29/combSys/comb_sys_edges_{var}_ch{region}.json"       ,"r") as f:
          comb_edges = json.load(f) 
   
        #build a th1d out of it
        comb_edges = array.array("d", comb_edges)
  
        hWeightsUp   = ROOT.TH1D("wHistUp","wHistUp", len(comb_up), comb_edges)
        hWeightsDown = ROOT.TH1D("wHistDown","wHistDown", len(comb_down), comb_edges)
    
        for i, w in enumerate(comb_up, start=1):
          hWeightsUp  .SetBinContent(i, w)
  
        for i, w in enumerate(comb_down, start=1):
          hWeightsDown.SetBinContent(i, w)
  
  
        ROOT.gInterpreter.Declare("""
        TH1D* hWUp = nullptr;

        double getWeightUp(double x) {
          return hWUp->GetBinContent(hWUp->FindBin(x));
        }
        """)
 
        ROOT.gInterpreter.Declare("""
        TH1D* hWDown = nullptr;

        double getWeightDown(double x) {
          return hWDown->GetBinContent(hWDown->FindBin(x));
        }
        """)
 
 
        ROOT.hWUp = hWeightsUp
        rdf_filt = rdf_filt.Define(f"comb_w_up"  , f"getWeightUp({var})")

        ROOT.hWDown = hWeightsDown
        rdf_filt = rdf_filt.Define(f"comb_w_down", f"getWeightDown({var})")

        print("====> Defined weights for combinatorial")


      ##############################################
      # Start filling                              #
      ##############################################

      if   var == "phiPi_m": tofill = "m_corr"
      elif var == "dxy_mu_err_pv": tofill = "dxy_mu_err_pv_corr"
      elif var == "dxy_mu_sig_pv": tofill = "dxy_mu_sig_pv_corr"
      elif var == "k1_pt": tofill = "k1_corr"
      elif var == "k2_pt": tofill = "k2_corr"
      else: tofill = var
      print(f"filling variable {tofill}")

      if (total_w_str == "" or data == True):

        histos[var] = rdf_filt\
                 .Histo1D(model[0], tofill)

      else:
        histos[var] = rdf_filt\
                 .Define("total_w", total_w_str)\
                 .Histo1D(model[0], tofill, "total_w")
        print(f"filling central curve with weight {total_w_str}")

      if ff_sys:


        if "star" not in mc: sys_dir = sys_scalar #only take e1 - e6 for scalar signals (BCL)
        else:                 sys_dir = sys_vector #take e1-e10

        for s in sys_dir:

          s_up   = s + "_up"
          s_down = s + "_down"

          var_up_av   = averages[s_up   + "_" + mc]
          var_down_av = averages[s_down + "_" + mc]

          total_w_s_up   = s_up
          total_w_s_down = s_down

          if mc:
            total_w_s_up   += " * trigger_sf"
            total_w_s_down += " * trigger_sf"


          ########################################
          # Note: Apply bs tau central here too? #
          ########################################

          if lt_central:

            total_w_s_up   += " * w_bs_tau"
            total_w_s_down += " * w_bs_tau"
 
          print(f"filling variational curve with weight {total_w_s_up}")
          print(f"filling variational curve with weight {total_w_s_down}")

 
          histos[var + "_" + s + "Up"]    = rdf_filt\
                                     .Define(f"total_w_{s_up}", total_w_s_up)\
                                     .Histo1D(model[0], tofill, f"total_w_{s_up}")
          histos[var + "_" + s + "Up"]    .Scale(1.0 / var_up_av)

          histos[var + "_" + s + "Down"]  = rdf_filt\
                                     .Define(f"total_w_{s_down}", total_w_s_down)\
                                     .Histo1D(model[0], tofill, f"total_w_{s_down}" )
          histos[var + "_" + s + "Down" ] .Scale(1.0 / var_down_av)


      if lt_sys:

          total_w_bs_up   = " w_bs_tau_up "
          total_w_bs_down = " w_bs_tau_down "

          if ff_central:
            total_w_bs_up   += f" * central_w  / {central_av}"
            total_w_bs_down += f" * central_w  / {central_av}"


          histos[var + "_bsTauUp"]    = rdf_filt\
                                          .Define(f"total_w_bs_up"    ,total_w_bs_up)\
                                          .Histo1D(model[0], tofill , "total_w_bs_up" )



          histos[var + "_bsTauDown"]  = rdf_filt\
                                          .Define(f"total_w_bs_down"    ,total_w_bs_down)\
                                          .Histo1D(model[0], tofill , "total_w_bs_down" )

      if combSys: 

          total_w_comb_up   = total_w_str 
          total_w_comb_down = total_w_str 

          total_w_comb_up   += f" * (comb_w_up  )"
          total_w_comb_down += f" * (comb_w_down)"

          print(total_w_comb_up)
      

          histos[var + "_combUp"]       = rdf_filt\
                                          .Define(f"total_w_comb_up"    ,total_w_comb_up)\
                                          .Histo1D(model[0], tofill , "total_w_comb_up" )



          histos[var + "_combDown"]     = rdf_filt\
                                          .Define(f"total_w_comb_down"    ,total_w_comb_down)\
                                          .Histo1D(model[0], tofill , "total_w_comb_down" )



      histos[var].GetXaxis().SetTitle(model[1])
      histos[var].SetMaximum(1.2 * histos[var].GetBinContent(histos[var].GetMaximumBin()))

      #get default colors 
      color,_ = getColorAndLabel(var)
      histos[var].SetLineColor(color)
      histos[var].SetLineWidth(2)


  return histos



def saveHisto1D(hist_dict, name):

  #extract root Th1d histo
  try   : histo1d_val = {name: h.GetValue() for name, h in hist_dict.items()}
  except: histo1d_val = {name: h            for name, h in hist_dict.items()}

  #save into file  
  f = ROOT.TFile(name, "RECREATE")
  for name, h in histo1d_val.items():
    h.SetName(name)   
    h.Write()
  f.Close()


#load fitted_vars
with open(f"{args.toSave_plots}/fitted_vars.json", "r") as f:
  fitted_vars = json.load(f)
  print(f" ====> JSON file loaded: {fitted_vars}")

with open(f"{args.toSave_plots}/normalization_parameters.json", "r") as f:
  norm_params = json.load(f)
  print(f" ====> JSON file loaded: {norm_params}")

global tau, n0_exp 
tau    = norm_params["tau"]
n0_exp = norm_params["n0_exp"]
m0     = norm_params["m0"]
m1     = norm_params["m1"]
tilt   = norm_params["tilt"]

#make this a class 
selec = selections(args.selection)

if args.channel == "DsMu":
  histo1d       = createHistos(selec.dsMu + hammer_str        ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsmu",      region = args.bin, lt_central = bs_tau_central, lt_sys = bs_tau_sys)
  print(histo1d)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_{args.bin}.root")

if args.channel == "DsMu_woHammer":
  histo1d       = createHistos(selec.dsMu + hammer_str        ,    rdfSig        ,                                                     mc = "dsmu",      region = args.bin, lt_central = bs_tau_central, lt_sys = bs_tau_sys) 
  print(histo1d)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_woHammer_{args.bin}.root")

if args.channel == "DsTau":

  histo1d       = createHistos(selec.dsTau + hammer_str       ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dstau",     region = args.bin, lt_central = bs_tau_central, lt_sys = bs_tau_sys)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsTau_{args.bin}.root")
  
  histo1d_blind = { key: histo1d[key].Clone()  for key in histo1d.keys() }
  for key in histo1d_blind.keys()    : histo1d_blind[key]    .Scale(blind_scalar)    
  saveHisto1D(histo1d_blind, f"{args.toSave_plots}/histos_DsTau_blind_{args.bin}.root")

if args.channel == "DsStarMu":
  histo1d       = createHistos(selec.dsStarMu + hammer_str    ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsstarmu",  region = args.bin, lt_central = bs_tau_central, lt_sys = bs_tau_sys)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsStarMu_{args.bin}.root")

if args.channel == "DsStarTau":
  histo1d       = createHistos(selec.dsStarTau + hammer_str   ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsstartau", region = args.bin, lt_central = bs_tau_central, lt_sys = bs_tau_sys)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsStarTau_{args.bin}.root")

  histo1d_blind = { key: histo1d[key].Clone()  for key in histo1d.keys() }
  for key in histo1d_blind.keys()    : histo1d_blind[key]    .Scale(blind_vector)    
  saveHisto1D(histo1d_blind, f"{args.toSave_plots}/histos_DsStarTau_blind_{args.bin}.root")

if args.channel == "Hb":
  histo1d       = createHistos(selec.hb           ,    rdfHb         ,                            mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_{args.bin}.root")

if args.channel == "Hb_dc":
  histo1d       = createHistos(selec.hb_dc        ,    rdfHb         ,                            mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_dc_{args.bin}.root")

if args.channel == "Hb_fd":
  histo1d       = createHistos(selec.hb_fd        ,    rdfHb         ,                            mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_fd_{args.bin}.root")



if args.channel == "Hb_others":
  histo1d       = createHistos(selec.hb_others        ,    rdfHb         ,                        mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_others_{args.bin}.root")



if args.channel == "Hb_bpm":
  histo1d       = createHistos(selec.hb_bpm        ,    rdfHb         ,                           mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bpm_{args.bin}.root")

if args.channel == "Hb_bpm_fd":
  histo1d       = createHistos(selec.hb_bpm_fd        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bpm_fd_{args.bin}.root")

if args.channel == "Hb_bpm_dc":
  histo1d       = createHistos(selec.hb_bpm_dc        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bpm_dc_{args.bin}.root")



if args.channel == "Hb_b0":
  histo1d       = createHistos(selec.hb_b0        ,    rdfHb         ,                            mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_b0_{args.bin}.root")

if args.channel == "Hb_b0_fd":
  histo1d       = createHistos(selec.hb_b0_fd        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_b0_fd_{args.bin}.root")

if args.channel == "Hb_b0_dc":
  histo1d       = createHistos(selec.hb_b0_dc        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_b0_dc_{args.bin}.root")



if args.channel == "Hb_bs":
  histo1d       = createHistos(selec.hb_bs        ,    rdfHb         ,                            mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bs_{args.bin}.root")

if args.channel == "Hb_bs_fd":
  histo1d       = createHistos(selec.hb_bs_fd        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bs_fd_{args.bin}.root")

if args.channel == "Hb_bs_dc":
  histo1d       = createHistos(selec.hb_bs_dc        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_bs_dc_{args.bin}.root")



if args.channel == "Hb_lambdab":
  histo1d       = createHistos(selec.hb_lambdab        ,    rdfHb         ,                       mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_lambdab_{args.bin}.root")

if args.channel == "Hb_lambdab_fd":
  histo1d       = createHistos(selec.hb_lambdab_fd        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_lambdab_fd_{args.bin}.root")

if args.channel == "Hb_lambdab_dc":
  histo1d       = createHistos(selec.hb_lambdab_dc        ,    rdfHb         ,                         mc = "hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_lambdab_dc_{args.bin}.root")



if args.channel == "Mu_in_Hb":
  histo1d       = createHistos(selec.dsMu              ,    rdfHb         ,                mc = "mu_in_hb", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Mu_in_Hb_{args.bin}.root")

if args.channel == "Data":
  histo1d       = createHistos(selec.bare         ,    rdfData       ,                                                   region = args.bin, data = True)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Data_{args.bin}.root")

if args.channel == "Data_sf_pimu":
  histo1d       = createHistos(args.data_selection,  rdfData         ,  sf_weights = bdt, sf_weights2 = bdt2,            region = args.bin, comb = "comb", combSys = True)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Data_sf_pimu_{args.bin}.root")


