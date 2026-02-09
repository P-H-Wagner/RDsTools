import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
#ROOT.ROOT.EnableImplicitMT() 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from sidebands import getSigma, getABCS
from signflip  import getSignflipRatio, getSignflipRatioTest, fitAnotherVar
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
parser.add_argument("--nn",          required = True,     help = "Specify 'before' or 'past' neural network plots ") 
parser.add_argument("--hammer",      required = True,     help = "Specify 'true' or 'false' to apply hammer weights") 
parser.add_argument("--hammer_sys",  required = True,     help = "Specify 'true' or 'false' to save weight variation shapes") 
parser.add_argument("--trigger",     required = True,     help = "Specify mu7 or mu9 for the trigger") 
parser.add_argument("--bdt",         required = True,     help = "Specify 'true' or 'false' to add bdt weights") 
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

print(f"====> Running  {args.nn} neural network on trigger {args.trigger} and sf_weights: {args.bdt} with hammer on: {args.hammer} and sys {args.hammer_sys}")

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


print(f"====> Running  {args.nn} neural network on trigger {args.trigger} and sf_weights: {args.bdt} with hammer on: {args.hammer} and sys {args.hammer_sys}")

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

with open( f"{args.toSave_plots}/info.txt", "a") as f:
  print(f"====> Running  {args.nn} neural network on trigger {args.trigger} and sf_weights: {args.bdt} with hammer on: {args.hammer} and sys {args.hammer_sys}")

#############
# Split MC  #
#############

if trigger == "mu7":
  data_selec   = " && (mu7_ip4 == 1)"
  mc_selec     = " && (mu7_ip4 == 1) && (static_cast<int>(event) % 20 < 10) " #on mc we use event nr (all triggers are on for mc!)
else:
  data_selec   = " && ((mu9_ip6 == 1) && (mu7_ip4 == 0)) "
  mc_selec     = " && (mu9_ip6 == 1) && (static_cast<int>(event) % 20 >= 10) "


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
  print("constrained data after NN not processed yet..")

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
      if var == "score1" and (args.split == "class") and region != None:
        
        #adapt the binning
        print(f"Adapt binning for {var} and region {region}")
  
        model = special_models[var + f"_bin{region}" ]

      #if var == "score1" and (args.split == "q2_coll" or args.split == "q2_lhcb_alt") and region != None:
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


def saveHisto1D(hist_dict, name):

  #extract root Th1d histo
  histo1d_val = {name: h.GetValue() for name, h in hist_dict.items()}

  #save into file  
  f = ROOT.TFile(name, "RECREATE")
  for name, h in histo1d_val.items():
    h.SetName(name)   
    h.Write()
  f.Close()



#make this a class 
selec = selections(args.selection)

if args.channel == "DsMu":
  histo1d       = createHistos(selec.dsMu         ,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsmu",      region = args.bin)
  print("created dsmu histos, saving now ....")
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_{args.bin}.root")

if args.channel == "DsMu_woHammer":
  histo1d       = createHistos(selec.dsMu         ,    rdfSig        , gen = False,                                                                               region = args.bin) 
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsMu_woHammer_{args.bin}.root")

if args.channel == "DsTau":
  histo1d       = createHistos(selec.dsTau        ,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dstau",     region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsTau_{args.bin}.root")

if args.channel == "DsStarMu":
  histo1d       = createHistos(selec.dsStarMu     ,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsstarmu",  region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsStarMu_{args.bin}.root")

if args.channel == "DsStarTau":
  histo1d       = createHistos(selec.dsStarTau    ,    rdfSig        , gen = False,  hammer_central = hammer_central, hammer_sys = hammer_sys, sig = "dsstartau", region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_DsStarTau_{args.bin}.root")

if args.channel == "Hb":
  histo1d       = createHistos(selec.hb           ,    rdfHb         , gen = False,                                                                                         region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Hb_{args.bin}.root")

if args.channel == "Mu_in_Hb":
  histo1d       = createHistos(selec.dsMu         ,    rdfHb         , gen = False,                                                                                         region = args.bin)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Mu_in_Hb_{args.bin}.root")

if args.channel == "Data":
  histo1d       = createHistos(selec.bare         ,    rdfData       , gen = False,                                                                                         region = args.bin, data = True)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Data_{args.bin}.root")

if args.channel == "Data_sf_pimu":
  histo1d       = createHistos(args.data_selection,  rdfData         , gen = False, sf_weights = sf_weights,                                                                region = args.bin, data = True)
  saveHisto1D(histo1d, f"{args.toSave_plots}/histos_Data_sf_pimu_{args.bin}.root")


