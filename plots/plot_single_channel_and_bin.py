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
from histModels import models, modelsSR, pastNN_models, pastNN_2Dmodels, special_models, special_models_q2_coll,special_models_e_star_lhcb_alt, special_models_ds_perp, special_models_score2
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
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
parser.add_argument("--bdt2",         required = True,     help = "Specify 'true' or 'false' to add bdt2 weights") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
parser.add_argument("--control",                          help = "If given, run control plots, either 'highmass', 'leftsb', 'rightsb' or 'complete' or 'custom' ") 
parser.add_argument("--cut",                              help = "Cut on the discriminator score 5. ") 
#parser.add_argument("--findcut",     action='store_true', help = "If given, we run thecut scan") 

parser.add_argument("--bin",   type=int, required = True,     help = "") 
parser.add_argument("--region",          required = True,     help = "") 
parser.add_argument("--toSave_plots",    required = True,     help = "") 
parser.add_argument("--selection",       required = True,     help = "") 
parser.add_argument("--data_selection",  required = True,     help = "") 
#parser.add_argument("--scale_pimu",      required = True,     help = "") 
#parser.add_argument("--scale_rest",      required = True,     help = "") 
parser.add_argument("--channel",         required = True,     help = "") 
parser.add_argument("--split",         required = True,     help = "") 

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

if args.control and args.control not in ["highmass","leftsb","rightsb", "complete", "custom"]:
    raise ValueError ("Error: Not a valid key for --control, please use 'highmass', 'sb' or 'complete' or 'custom' ")
control = args.control

if args.cut and not pastNN : raise ValueError ("Error: Cannot interpret cut before NN")
elif args.cut and pastNN   : score_cut = f" && (score5 <= {args.cut} ) "
#elif args.cut and pastNN   : score_cut = f" && (fv_prob > 0.1) && (rel_iso_03_pv < 0.3) && (bs_pt_lhcb_alt > 8) && (abs(cosPiK1) > 0.9) && (cosMuW_lhcb_alt< -0.2) && (bs_pt_lhcb_alt > 20) && (bs_mass_corr < 7) "
else                       : score_cut = ""

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
  f.write(f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} on nn {nn_model}\n")
  print  (f"====> Running  {args.nn} neural network and bdt: {bdt} and bdt2: {bdt2} with hammer on: {hammer_central} and sys {hammer_sys} on nn {nn_model}\n")



#############
# Split MC  #
#############

#if trigger == "mu7":
#  data_selec   = " && (mu7_ip4 == 1)"
#  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
#else:
#  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
#  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "

#data_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
data_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#data_selec = " && ((mu7_ip4 == 1)||(mu9_ip6==1))"
#data_selec = " && ((mu7_ip4 == 1))"
#data_selec = " && ((mu8_ip3 == 1))"
#data_selec = " && ((mu8_ip5 == 1))"
#data_selec = " && ((mu8_ip6 == 1))"
#data_selec = " && ((mu9_ip4 == 1))"
#data_selec = " && ((mu9_ip5 == 1))"
#data_selec = " && ((mu9_ip6 == 1))"
#data_selec = " && ((mu12_ip6 == 1))"


#mc_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip4==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = " && ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
mc_selec = "&& ((mu7_ip4 == 1)||(mu8_ip3==1)||(mu8_ip5==1)||(mu8_ip6==1)||(mu9_ip5==1)||(mu9_ip6==1)||(mu12_ip6))"
#mc_selec = " && ((mu7_ip4 == 1))"
#mc_selec = " && ((mu8_ip3 == 1))"
#mc_selec = " && ((mu8_ip5 == 1))"
#mc_selec = " && ((mu8_ip6 == 1))"
#mc_selec = " && ((mu9_ip4 == 1))"
#mc_selec = " && ((mu9_ip5 == 1))"
#mc_selec = " && ((mu9_ip6 == 1))"
#mc_selec  = " && ((mu12_ip6 == 1))"


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

#update data
if (bdt and not bdt2):

  files_data = [bdt_data]
  path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/"

#update sig and data (includes hammer weights and bdt weights anyways)
if pastNN:

  files_sig  = [sig_pastNN ]
  files_hb   = [hb_pastNN  ]
  #files_data = [data_pastNN]

  path_sig   = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees_test/"
  path_hb    = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees_test/"
  #path_data  = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/"

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


def createHistos(selection, rdf, data = False , variables = None, ff_central = False, ff_sys = False, mc = None, sf_weights = None, sf_weights2 = None, region = None, massfit = False):

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

  print(f"====> Total applied weight is: {total_w_str}" )

  for var,model in models.items():


      #in the case we have a special binning, only plot the corresponding variable for the region!
      if ((len(fitted_vars) != 0) and (var != fitted_vars[str(region)])): continue

      ################################
      # Special binnings             #
      ################################

      if var == "score1" and (args.split == "class") and region != None:

        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")

        model = special_models[var + f"_bin{region}" ]

      if var == "score1" and (args.split == "q2_coll" or args.split == "q2_lhcb_alt") and region != None:

        #adapt the binning
        print(f"======> Adapt binning for {var} and region {region}")

        model = special_models_q2_coll[var + f"_bin{region}" ]

      if var == "e_star_lhcb_alt" and (args.split == "q2_coll" or args.split == "q2_lhcb_alt") and region != None:

        #adapt the binning
        print(f"======> Adapt binning for {var} and region {region}")

        model = special_models_e_star_lhcb_alt[var + f"_bin{region}" ]


      if ("score" in var or var == "phiPi_m" or "q2" in var) and (args.split == "score2" ) and region != None:
        
        #adapt the binning
        print(f"======> Adapt binning for {var} and region {region}")
  
        model = special_models_score2[var + f"_bin{region}" ]  


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
      #mass_expr   = f"(run==1) * 1.000 * phiPi_m + (run!=1) * phiPi_m"

      if var == "phiPi_m": tofill = "m_corr"
      else: tofill = var

      if (total_w_str == "" or data == True):
        histos[var] = rdf.Filter(selection).Define("m_corr",mass_expr).Histo1D(model[0], tofill)
      else:
        histos[var] = rdf.Filter(selection).Define("m_corr",mass_expr).Define("total_w", total_w_str).Histo1D(model[0], tofill, "total_w")

      if ff_sys:

        #this is only true for signals
        #add also systematical shape variations as variables 

        if "star" not in mc: sys_dir = sys_scalar #only take e1 - e6 for scalar signals (BCL)
        else:                 sys_dir = sys_vector #take e1-e10

        for s in sys_dir:

          s_up   = s + "_up"
          s_down = s + "_down"

          var_up_av   = averages[s_up   + "_" + mc]
          var_down_av = averages[s_down + "_" + mc]

          func_up   = s_up   + f" / ({ var_up_av   })"
          func_down = s_down + f" / ({ var_down_av })"

          total_w_s_up   = s_up
          total_w_s_down = s_down
          if mc:
            total_w_s_up   += " * trigger_sf"
            total_w_s_down += " * trigger_sf"


          # fill histogram with weight "s"
 
          histos[var + "_" + s + "Up"]    = rdf.Filter(selection).Define("m_corr", mass_expr).Define(f"total_w_{s_up}", total_w_s_up).Histo1D(model[0], var, f"total_w_{s_up}")
          histos[var + "_" + s + "Up"]    .Scale(1.0 / var_up_av)

          histos[var + "_" + s + "Down"]  = rdf.Filter(selection).Define("m_corr", mass_expr).Define(f"total_w_{s_down}", total_w_s_down).Histo1D(model[0], var, f"total_w_{s_down}" )
          histos[var + "_" + s + "Down" ] .Scale(1.0 / var_down_av)


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
 

#make this a class 
selec = selections(args.selection)

if args.channel == "DsMu":
  histo1d       = createHistos(selec.dsMu + hammer_str        ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsmu",      region = args.bin)
  print(histo1d)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_{args.bin}.root")

if args.channel == "DsMu_woHammer":
  histo1d       = createHistos(selec.dsMu + hammer_str        ,    rdfSig        ,                                                     mc = "dsmu",      region = args.bin) 
  print(histo1d)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_woHammer_{args.bin}.root")

if args.channel == "DsTau":

  histo1d       = createHistos(selec.dsTau + hammer_str       ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dstau",     region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsTau_{args.bin}.root")
  
  histo1d_blind = { key: histo1d[key].Clone()  for key in histo1d.keys() }
  for key in histo1d_blind.keys()    : histo1d_blind[key]    .Scale(blind_scalar)    
  saveHisto1D(histo1d_blind, f"{args.toSave_plots}/histos_DsTau_blind_{args.bin}.root")

if args.channel == "DsStarMu":
  histo1d       = createHistos(selec.dsStarMu + hammer_str    ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsstarmu",  region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsStarMu_{args.bin}.root")

if args.channel == "DsStarTau":
  histo1d       = createHistos(selec.dsStarTau + hammer_str   ,    rdfSig        ,   ff_central = hammer_central, ff_sys = hammer_sys, mc = "dsstartau", region = args.bin)
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
  histo1d       = createHistos(args.data_selection,  rdfData         ,  sf_weights = bdt, sf_weights2 = bdt2,            region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Data_sf_pimu_{args.bin}.root")


