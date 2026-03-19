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
chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.datetime}/*")
df = ROOT.RDataFrame(chain)

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")


dest_dir = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dt}/"
os.system(f"mkdir -p {dest_dir}")

#calculate time t
ROOT.gInterpreter.Declare(r'''double get_time(float pv_x, float pv_y, float pv_z, float sv_x, float sv_y, float sv_z, float beta){

  //remark sv and pv quantities are in cm

  float distance = sqrt(pow(sv_x - pv_x,2) + pow(sv_y - pv_y,2) + pow(sv_z - pv_z,2));

  //convert distance to meters
  distance *= 0.01; 

  //remark boost is in atural units, i.e. a number in [0,1] (the beta)
  float gamma = 1.0 / sqrt(1 - pow(beta,2));

  // use c in SI units to get a time in seconds
  float c = 299792458.0;
  float t = distance / (gamma * beta * c);

  return t;
  }'''
  )

df.Define  ("time","get_time(gen_pv_x, gen_pv_y, gen_pv_z, gen_sv_x, gen_sv_y, gen_sv_z, gen_bs_boost)") \
  .Snapshot("tree",f"{dest_dir}/tree_with_time.root")
