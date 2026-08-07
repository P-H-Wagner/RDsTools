import ROOT
import numpy as np
import argparse
from datetime import datetime
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))


parser = argparse.ArgumentParser()
parser.add_argument("--datetime", required = True)
args = parser.parse_args()


#create RDF 
#input file is a gen lvl sample wout filter
chain = ROOT.TChain("tree")
chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.datetime}/*.root")
df = ROOT.RDataFrame(chain)

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")

dest_dir = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dt}/"
os.system(f"mkdir -p {dest_dir}")

#calculate time t
ROOT.gInterpreter.Declare(r'''double get_rapidity(float bs_eta, float bs_phi, float bs_pt, float bs_m){

  //change coords to x y z
  float p_x = bs_pt * cos (bs_phi);
  float p_y = bs_pt * sin (bs_phi);
  float p_z = bs_pt * sinh(bs_eta);

  float e   = sqrt(pow(p_x,2) + pow(p_y,2) + pow(p_z,2) + pow(bs_m,2));

  float y   = 0.5 * log((e + p_z) / (e - p_z));

  return y;
  }'''
  )

df.Define  ("bs_y","get_rapidity(gen_bs_eta, gen_bs_phi, gen_bs_phi, gen_bs_m)") \
  .Snapshot("tree",f"{dest_dir}/tree_with_y.root")
