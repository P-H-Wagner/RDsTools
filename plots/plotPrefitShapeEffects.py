import ROOT
import argparse
import os
import sys
import yaml
from datetime import datetime
import matplotlib.pyplot as plt
import glob
import pdb

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/comb"))
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))

ROOT.gROOT.SetBatch(True)

# parsing
parser = argparse.ArgumentParser()
parser.add_argument("--file",         required = True,     help = "Specify dt of file")
args = parser.parse_args()

include_lifetime=False

#inside here we have all the root files, containing all signal shapes + up/down variation
shapes_dir= f"/work/pahwagne/RDsTools/fit/shapes_binned/{args.file}/*"
files = glob.glob(shapes_dir) 

#save at
dest = f"/work/pahwagne/RDsTools/plots/prefitShapeVariations/{args.file}"
os.system(f"mkdir -p {dest}")


var     = set() 
regions = set()   

for f in files:

  #pdb.set_trace()
  #pick all root files
  if ".root" in f:
    #remove .root
    f = f.split(".root")[0]
    #split 
    #var_bin  = f.split("_shapes_constrained_pastNN_")
    var_bin  = f.split("_shapes_")

    var    .add(var_bin[0])
    regions.add(var_bin[1])

    #print(var_bin)


def createCanvas(name):


  c = ROOT.TCanvas(name,name, 800, 800)

  # divide canvas for different signals
  #if ((name == "dsmu") or (name == "dstau")):
  #  c = ROOT.TCanvas(name,name, 20000, 30000)
  #  #bcl has 6
  #  ncols = 2
  #  nrows = 3

  #  c.Divide(ncols, nrows)

  #if ((name == "dsstmu") or (name == "dssttau")):
  #  c = ROOT.TCanvas(name,name, 20000, 60000)
  #  #bgl has 10
  #  ncols = 2
  #  nrows = 5

  #  c.Divide(ncols, nrows)

  return c

def plotHammer(rf, sig, var, chan):

  if "Star" in sig: 
    model = "Bgl"
    nParams = 11
  else: 
    model = "Bcl"
    nParams = 7

  os.system(f"mkdir -p {dest}/{sig}")
  for i in range(1,nParams):

    c = createCanvas(f"{sig}_e{i}")
    c.Draw()
    c.cd()

    main_pad = ROOT.TPad('main_pad', '', 0., 0.25, 1. , 1.  )
    main_pad.Draw()
    ratio_pad = ROOT.TPad('ratio_pad', '', 0., 0., 1., 0.25)
    ratio_pad.Draw()

    main_pad.SetTicks(True)
    main_pad.SetBottomMargin(0.)
    main_pad.SetLeftMargin(.16)

    ratio_pad.SetTopMargin(0.)
    ratio_pad.SetLeftMargin(.16)
    ratio_pad.SetGridy()
    ratio_pad.SetBottomMargin(0.45)

    #############################
    main_pad.cd()

    hUp   = rf.Get(f"{sig}_e{i}{model}Up_ch{chan}")
    hUp.SetLineColor(ROOT.kRed)
    hC    = rf.Get(f"{sig}_ch{chan}")
    hC.SetLineColor(ROOT.kBlack)
    hDown = rf.Get(f"{sig}_e{i}{model}Down_ch{chan}")
    hDown.SetLineColor(ROOT.kBlue)
    
    hUp.SetTitle(f"{sig}")
    hUp.GetXaxis().SetTitle(var)
    hUp.GetYaxis().SetTitle("entries")

    hUp   .SetMinimum(0)
    hUp.GetYaxis().SetRangeUser(1e-3, hUp.GetBinContent(hUp.GetMaximumBin())*1.5)
    hUp  .Draw("E")    
    hC   .Draw("E SAME")    
    hDown.Draw("E SAME")    


    legend = ROOT.TLegend(0.6,0.7,0.9,0.85)
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hUp  , fr"{model} e{i}Up in ch{chan}", "L")  
    legend.AddEntry(hDown, fr"{model} e{i}Down in ch{chan}", "L")  
    legend.AddEntry(hC   , f"Central", "L")  
    legend.Draw("SAME")

  
    ############################# 
    ratio_pad.cd()

    ratio_up = hUp.Clone()
    ratio_up.Divide(hC)
    ratio_up.SetLineColor(ROOT.kRed)
    ratio_up.SetTitle("")
    ratio_up.GetYaxis().SetTitle("Up(Down) / Central")

    ratio_down = hDown.Clone()
    ratio_down.Divide(hC)
    ratio_down.SetLineColor(ROOT.kBlue)

    ratio_up.GetYaxis().SetNdivisions(505)
    ratio_up.GetYaxis().SetTitleSize(0.1)
    ratio_up.GetYaxis().SetLabelSize(0.08)
    ratio_up.GetXaxis().SetTitleSize(0.1)
    ratio_up.GetXaxis().SetLabelSize(0.08)

    ratio_max = max([ratio_up.GetMaximum(), ratio_down.GetMaximum()])
    ratio_min = min([ratio_up.GetMinimum(), ratio_down.GetMinimum()])

    ratio_up.SetMinimum(0.9)
    ratio_up.SetMaximum(1.099)
    ratio_up.GetYaxis().SetTitleOffset(0.5)

    ROOT.gPad.RedrawAxis()
    
    ratio_up.Draw("E")
    ratio_down.Draw("E SAME")
 
    c.Modified()
    c.Update()
    c    .SaveAs(f"{dest}/{sig}/{var}_in_bin_ch{chan}_{sig}_e{i}_variations.pdf")

    ## check if it is a true envelope
    nbins = hUp.GetNbinsX()
    for j in range(1,nbins+1):
      nUp   = int(hUp  .GetBinContent(j))
      nDown = int(hDown.GetBinContent(j))
      nC    = int(hC   .GetBinContent(j))

      if ((nC < nDown) and (nC < nUp)) or ((nC > nDown) and (nC > nUp)): 
        print(f"ALERT !!! - There is some one-sided variation in signal {sig} and channel e{i} - abort")
        print(f"Bin content: \nup = {nUp} \ndown = {nDown} \ncentral = {nC}")
        sys.exit()


def plotLifetime(rf,sig, var):

    #lifetime
    c = createCanvas(f"{sig}_bs_tau")
    c.cd()
    hUp   = rf.Get(f"{sig}_bsTauUp_ch{chan}")
    hUp.SetLineColor(ROOT.kRed)
    hC    = rf.Get(f"{sig}_ch{chan}")
    hC.SetLineColor(ROOT.kBlack)
    hDown = rf.Get(f"{sig}_bsTauDown_ch{chan}")
    hDown.SetLineColor(ROOT.kBlue)
  
    hUp.SetTitle(f"{sig}")
    hUp.GetXaxis().SetTitle(var)
    hUp.GetYaxis().SetTitle("entries")
  
    hUp   .SetMinimum(0)
    hUp.GetYaxis().SetRangeUser(1e-3, hUp.GetBinContent(hUp.GetMaximumBin())*1.5)
    hUp  .Draw("E")    
    hC   .Draw("E SAME")    
    hDown.Draw("E SAME")    
   
    legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hUp  , f"bsTauUp in ch{chan}", "L")  
    legend.AddEntry(hDown, f"bsTauDown in ch{chan}", "L")  
    legend.AddEntry(hC   , f"Central", "L")  
    legend.Draw("SAME")
  
    c    .SaveAs(f"{dest}/{sig}/{var}_in_bin_ch{chan}_{sig}_bsTau_variations.pdf")

def plotCombSys(rf,sig, var, chan):

    print(chan)
    os.system(f"mkdir -p {dest}/{sig}")

    #lifetime
    c = createCanvas(f"{sig}_comb")
    c.cd()
    hUp   = rf.Get(f"{sig}_combSysUp_ch{chan}")
    hUp.SetLineColor(ROOT.kRed)
    hC    = rf.Get(f"{sig}_ch{chan}")
    hC.SetLineColor(ROOT.kBlack)
    hDown = rf.Get(f"{sig}_combSysDown_ch{chan}")
    hDown.SetLineColor(ROOT.kBlue)
  
    hUp.SetTitle(f"{sig}")
    hUp.GetXaxis().SetTitle(var)
    hUp.GetYaxis().SetTitle("entries")
  
    hUp   .SetMinimum(0)
    hUp.GetYaxis().SetRangeUser(1e-3, hUp.GetBinContent(hUp.GetMaximumBin())*1.5)
    hUp  .Draw("E")    
    hC   .Draw("E SAME")    
    hDown.Draw("E SAME")    
   
    legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
    legend.SetTextSize(0.04)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hUp  , f"combUp in ch{chan}", "L")  
    legend.AddEntry(hDown, f"combDown in ch{chan}", "L")  
    legend.AddEntry(hC   , f"Central", "L")  
    legend.Draw("SAME")
  
    c    .SaveAs(f"{dest}/{sig}/{var}_in_bin_ch{chan}_{sig}_comb_variations.pdf")



#for every bin/region, we have 4 canvas', one for each signal, showing all 6(10) variations

import pdb 
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
  var       = f.split("_shapes_")[0]
  split_ch  = f.split("_shapes_")[1]
  chan      = split_ch.split("ch")[1]
 

  for k in keys:
    key = k.GetName() 


  plotHammer(rf,"dsMu",      var, chan)
  plotHammer(rf,"dsStarMu",  var, chan)
  plotHammer(rf,"dsTau",     var, chan)
  plotHammer(rf,"dsStarTau", var, chan)

  plotCombSys  (rf,"comb", var, chan)


  if include_lifetime == True:

  
    plotLifetime(rf,"dsMu",      var)
    plotLifetime(rf,"dsStarMu",  var)
    plotLifetime(rf,"dsTau",     var)
    plotLifetime(rf,"dsStarTau", var)


