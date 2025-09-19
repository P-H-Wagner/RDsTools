import ROOT
from ROOT import TLorentzVector

import os
import argparse
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import pdb
ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!


sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

from helper import * 
import numpy as np

#ROOT.ROOT.EnableImplicitMT(8)

def boolean_string(s):
    if s not in {"false", "true"}:
        raise ValueError("Error: Not a valid boolean string, please use 'true' or 'false' (all lowercase!)")
    return s == "true"

# parsing
parser = argparse.ArgumentParser()
#parser.add_argument("--prod",        required = True,     help = "Specify '24' or '25' to specify the data production") 
#parser.add_argument("--trigger",     required = True,     help = "Specify mu7 or mu9 for the trigger") 
parser.add_argument("--debug",       action='store_true', help = "If given, run plotter with 50k events only") 
args = parser.parse_args()


######################################

#load the data

chain = ROOT.TChain("tree")
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/"
if args.debug:
  files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/bdt_weighted_data/{bdt_data_25}/*22*"

chain.Add(files)

# turn into rdf
rdf = ROOT.RDataFrame(chain)

# one permill scale
scale = 0.999

ROOT.gROOT.ProcessLine(r'''TLorentzVector scaleP4(double scale, double pt, double eta, double phi, double m){

  TLorentzVector p;
  p.SetPtEtaPhiM(scale*pt,eta,phi,m);
  return p;
  }
''')


ROOT.gROOT.ProcessLine(r'''double get_q2(TLorentzVector a, TLorentzVector b){

  TLorentzVector q = a - b;
  return q.M2();
  }
''')



ROOT.gROOT.ProcessLine(r'''TLorentzVector get_bs_coll(TLorentzVector ds, TLorentzVector mu){

  TLorentzVector bs = ds + mu;

  return bs *= 5.36688 / bs.M();
 
  }
''')


df = rdf.Define("ds_p4_scaled"     , "scaleP4(0.999, ds_refitted_pt     , ds_refitted_eta     , ds_refitted_phi  , 1.96834)") \
        .Define("mu_p4_scaled"     , "scaleP4(0.999, mu_refitted_pt     , mu_refitted_eta     , mu_refitted_phi  , 0.105658)") \
        .Define("bs_coll_p4_scaled", "get_bs_coll(ds_p4_scaled, mu_p4_scaled)") \
        .Define("q2_scaled"        , "get_q2(ds_p4_scaled, bs_coll_p4_scaled)") \
        #.Define("q2_unfitted"      , "get_q2(ds_p4_unfitted, bs_coll_p4_scaled)")
        #.Define("ds_p4_unfitted"   , "scaleP4(0.999, phiPi_pt           , phiPi_eta           , phiPi_phi        , 1.96834)") \


h = df.Histo1D(("q2_scaled", "", 31, 0, 12), "q2_scaled")
h.SetLineColor(ROOT.kBlue) 
h.GetYaxis().SetRangeUser(0,h.GetBinContent(h.GetMaximumBin()) * 1.5)
h2 = df.Histo1D(("q2_coll" , "", 31, 0, 12), "q2_coll") 
h.SetLineColor(ROOT.kRed) 

legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
legend.SetTextSize(0.04)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(h.GetValue(), "q2 permill scaled", "L")
legend.AddEntry(h2.GetValue(), "q2 unscaled", "L")


c = ROOT.TCanvas("c","c",800,600)
h.Draw()
h2.Draw("SAME")
legend.Draw("SAME")
c.Draw()
c.SaveAs("q2_permill_shift.png")
