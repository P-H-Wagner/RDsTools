import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import glob

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))


# parsing
parser = argparse.ArgumentParser()
parser.add_argument("--file",         required = True,     help = "Specify dt of file")
args = parser.parse_args()


#inside here we have all the root files, containing all signal shapes + up/down variation
shapes_dir= f"/work/pahwagne/RDsTools/fit/shapes_binned/{args.file}/*"
files = glob.glob(shapes_dir) 

#save at
dest = f"/work/pahwagne/RDsTools/plots/prefitShapeVariations/{args.file}"
os.system(f"mkdir -p {dest}")


var     = set() 
regions = set()   

for f in files:

  #pick all root files
  if ".root" in f:
    #remove .root
    f = f.split(".root")[0]
    #split 
    var_bin  = f.split("_shapes_constrained_pastNN_")

    var    .add(var_bin[0])
    regions.add(var_bin[1])

    print(var_bin)


def createCanvas(name):

  width  = 1200 
  height = 800

  c = ROOT.TCanvas(name,name, width, height)

  # divide canvas for different signals
  if ((name == "dsmu") or (name == "dstau")):
    #bcl has 6
    ncols = 2
    nrows = 3

    c.Divide(ncols, nrows)

  if ((name == "dsstmu") or (name == "dssttau")):
    #bgl has 10
    ncols = 2
    nrows = 5

    c.Divide(ncols, nrows)

  return c


#for every bin/region, we have 4 canvas', one for each signal, showing all 6(10) variations

 
for f in files:

  if ".root" not in f: continue;

  print(f"Producing plots for file {f}")
  rf = ROOT.TFile.Open(f, "READ") 
  keys = rf.GetListOfKeys()

  #check if list is not emptyy (can happen f.e. for 0 event bins!)
  if len(keys) <= 0: continue;

  #drop .root
  f = f.split(".root")[0]
  #remove path
  f = f.split(f"{args.file}/")[1]
  #split into variable and bin 
  var_bin   = f.split("_shapes_constrained_pastNN_")
 
  #create canvas
  c_0  = createCanvas("dsmu"   )
  c_1  = createCanvas("dstau"  )
  c_10 = createCanvas("dsstmu" )
  c_11 = createCanvas("dssttau")
 
  for k in keys:
   
    key = k.GetName() 
    #first, get central curves and draw them in every pad
    if (("dsMu"      in key) and ("Up" not in key) and ("Down" not in key)): central_0  = key
    if (("dsTau"     in key) and ("Up" not in key) and ("Down" not in key)): central_1  = key
    if (("dsStarMu"  in key) and ("Up" not in key) and ("Down" not in key)): central_10 = key
    if (("dsStarTau" in key) and ("Up" not in key) and ("Down" not in key)): central_11 = key

  for i in range(6):

    #draw central, scalar shapes
    c_0.cd(i+1)
    h = rf.Get(central_0)
    h.SetFillStyle(0)
    h.SetLineWidth(1)

    text = ROOT.TPaveText(.18,.7,.5,.88,"NDC");
    text.SetFillColor(0)
    text.SetFillStyle(0)
    text.SetTextAlign(11)
    text.AddText(f"Direction e{i+1}")
    h.Draw("HIST")
    text.Draw("SAME")

  #  c_1.cd(i+1)
  #  h = rf.Get(central_1)
  #  h.SetFillStyle(0)
  #  h.SetLineWidth(1)
  #  h.Draw("HIST")

  #for i in range(10):

  #  #draw central, vector shapes
  #  c_10.cd(i+1)
  #  h = rf.Get(central_10)
  #  h.SetFillStyle(0)
  #  h.SetLineWidth(1)
  #  h.Draw("HIST")

  #  c_11.cd(i+1)
  #  h = rf.Get(central_11)
  #  h.SetFillStyle(0)
  #  h.SetLineWidth(1)
  #  h.Draw("HIST")

  for k in keys:
   
    key = k.GetName() 

    if ("Up" not in key) and ("Down" not in key)          : continue;
    if ("hb" in key) or ("data" in key) or ("comb" in key): continue;

    print("Drawing variation for key: ", key)

    #draw variational shape
    signal = key.split("_e")[0] #f.e. dsStarMu
    rest   = key.split("_e")[1] #f.e. 8BglUp_ch0
    direc  = rest.split("B")[0] 

    if signal == "dsMu"     :
      print("This is dsmu in pad: ", int(direc)+1)
      c_0.cd(int(direc)+1)
 
  #  if signal == "dsTau"    : 
  #    print("This is dstau in pad: ", int(direc)+1)
  #    c_1.cd(int(direc)+1)

  #  if signal == "dsStarMu" : 
  #    print("This is ds*mu in pad: ", int(direc)+1)
  #    c_10.cd(int(direc)+1)

  #  if signal == "dsStarTau": 
  #    print("This is ds*tau in pad: ", int(direc)+1)
  #    c_11.cd(int(direc)+1)

      h = rf.Get(key)
      h.SetLineWidth(1)
      h.Draw("HIST SAME")

  c_0 .SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsmu_variations.pdf")
  c_1 .SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dstau_variations.pdf")
  c_10.SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsstarmu_variations.pdf")
  c_11.SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsstartau_variations.pdf")

 

  
