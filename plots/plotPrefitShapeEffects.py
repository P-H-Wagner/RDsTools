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

ROOT.gROOT.SetBatch(True)

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



  # divide canvas for different signals
  if ((name == "dsmu") or (name == "dstau")):
    c = ROOT.TCanvas(name,name, 20000, 30000)
    #bcl has 6
    ncols = 2
    nrows = 3

    c.Divide(ncols, nrows)

  if ((name == "dsstmu") or (name == "dssttau")):
    c = ROOT.TCanvas(name,name, 20000, 60000)
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


  titles = []
  for i in range(6):

    c_0.cd()
    #draw title on main canvas
    title = ROOT.TPaveText(0.3, 0.94, 0.7, 0.99, "NDC")  # top center
    title.AddText(r"B_{s} #rightarrow D_{s} #mu #nu")
    title.SetTextAlign(22)
    title.SetTextSize(0.02)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    titles.append(title)
    title.Draw("SAME")  

    #draw central, scalar shapes
    c_0.cd(i+1)
    h = rf.Get(central_0)
    h.SetName(f"dsMu_central_bin{i+1}")
    #h.SetFillStyle(0)
    h.SetLineWidth(1)
    h.Draw("E")
    h.GetYaxis().SetRangeUser(1e-3, h.GetBinContent(h.GetMaximumBin())*1.2)
    c_0.Update() 
 
    c_1.cd()
    #draw title on main canvas
    title = ROOT.TPaveText(0.3, 0.94, 0.7, 0.99, "NDC")  # top center
    title.AddText(r"B_{s} #rightarrow D_{s} #tau #nu")
    title.SetTextAlign(22)
    title.SetTextSize(0.02)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    titles.append(title)
    title.Draw("SAME")  

    #draw central, scalar shapes
    c_1.cd(i+1)
    h = rf.Get(central_1)
    h.SetName(f"dsMu_central_bin{i+1}")
    #h.SetFillStyle(0)
    h.SetLineWidth(1)
    h.Draw("E")
    h.GetYaxis().SetRangeUser(1e-3, h.GetBinContent(h.GetMaximumBin())*1.2)
    c_1.Update() 

  for i in range(10): 
    c_10.cd()
    #draw title on main canvas
    title = ROOT.TPaveText(0.3, 0.94, 0.7, 0.99, "NDC")  # top center
    title.AddText(r"B_{s} #rightarrow D*_{s} #mu #nu")
    title.SetTextAlign(22)
    title.SetTextSize(0.02)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    titles.append(title)
    title.Draw("SAME")  

    #draw central, scalar shapes
    c_10.cd(i+1)
    h = rf.Get(central_10)
    h.SetName(f"dsMu_central_bin{i+1}")
    #h.SetFillStyle(0)
    h.SetLineWidth(1)
    h.Draw("E")
    h.GetYaxis().SetRangeUser(1e-3, h.GetBinContent(h.GetMaximumBin())*1.2)
    c_10.Update()   

    c_11.cd()
    #draw title on main canvas
    title = ROOT.TPaveText(0.3, 0.94, 0.7, 0.99, "NDC")  # top center
    title.AddText(r"B_{s} #rightarrow D*_{s} #tau #nu")
    title.SetTextAlign(22)
    title.SetTextSize(0.02)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    titles.append(title)
    title.Draw("SAME")  

    #draw central, scalar shapes
    c_11.cd(i+1)
    h = rf.Get(central_11)
    h.SetName(f"dsMu_central_bin{i+1}")
    #h.SetFillStyle(0)
    h.SetLineWidth(1)
    h.Draw("E")
    h.GetYaxis().SetRangeUser(1e-3, h.GetBinContent(h.GetMaximumBin())*1.2)
    c_11.Update()   


  legends = []
  for k in keys:

    #ROOT.gPad.SetLogy()   
    key = k.GetName() 

    if ("Up" not in key) and ("Down" not in key)          : continue;
    if ("hb" in key) or ("data" in key) or ("comb" in key): continue;

    print("Drawing variation for key: ", key)

    #draw variational shape
    signal = key.split("_e")[0] #f.e. dsStarMu
    rest   = key.split("_e")[1] #f.e. 8BglUp_ch0
    direc  = rest.split("B")[0] 

    if signal == "dsMu"     :
      print("This is dsmu in pad: ", int(direc)) #no need for +1 here, we extract the nr from the name (starting at 1 already)
      c_0.cd(int(direc))
      c_0.Update() 

      import pdb
      #pdb.set_trace() 
      h = rf.Get(key)
      h.SetLineWidth(1)
      #h.SetLineStyle(2)
      h.SetLineColor(ROOT.kRed)
      h.Draw("E SAME")
      c_0.Update() 

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)
      #extract central curve
      pad = ROOT.gPad 
      objs = pad.GetListOfPrimitives()
      for obj in objs:
        print(obj.GetName(), "->", obj.ClassName())
        if "central" in obj.GetName(): 
          legend.AddEntry(obj, "Central curve", "L")  
      legend.AddEntry(h, f"Deviation e{direc}", "L")  
      legend.Draw("SAME")
      legends.append(legend) 
      c_0.Update() 


    if signal == "dsTau"    : 
      print("This is dstau in pad: ", int(direc))
      c_1.cd(int(direc))
      import pdb
      #pdb.set_trace() 
      h = rf.Get(key)
      h.SetLineWidth(1)
      #h.SetLineStyle(2)
      h.SetLineColor(ROOT.kRed)
      h.Draw("E SAME")
      c_1.Update() 

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)
      #extract central curve
      pad = ROOT.gPad 
      objs = pad.GetListOfPrimitives()
      for obj in objs:
        print(obj.GetName(), "->", obj.ClassName())
        if "central" in obj.GetName(): 
          legend.AddEntry(obj, "Central curve", "L")  
      legend.AddEntry(h, f"Deviation e{direc}", "L")  
      legend.Draw("SAME")
      legends.append(legend) 
      c_1.Update() 

    if signal == "dsStarMu" : 
      print("This is ds*mu in pad: ", int(direc))
      c_10.cd(int(direc))

      h = rf.Get(key)
      h.SetLineWidth(1)
      #h.SetLineStyle(2)
      h.SetLineColor(ROOT.kRed)
      h.Draw("E SAME")
      c_10.Update() 

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)
      #extract central curve
      pad = ROOT.gPad 
      objs = pad.GetListOfPrimitives()
      for obj in objs:
        print(obj.GetName(), "->", obj.ClassName())
        if "central" in obj.GetName(): 
          legend.AddEntry(obj, "Central curve", "L")  
      legend.AddEntry(h, f"Deviation e{direc}", "L")  
      legend.Draw("SAME")
      legends.append(legend) 
      c_10.Update() 

    if signal == "dsStarTau": 
      print("This is ds*tau in pad: ", int(direc))
      c_11.cd(int(direc))

      #pdb.set_trace() 
      h = rf.Get(key)
      h.SetLineWidth(1)
      #h.SetLineStyle(2)
      h.SetLineColor(ROOT.kRed)
      h.Draw("E SAME")
      c_11.Update() 

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)
      #extract central curve
      pad = ROOT.gPad 
      objs = pad.GetListOfPrimitives()
      for obj in objs:
        print(obj.GetName(), "->", obj.ClassName())
        if "central" in obj.GetName(): 
          legend.AddEntry(obj, "Central curve", "L")  
      legend.AddEntry(h, f"Deviation e{direc}", "L")  
      legend.Draw("SAME")
      legends.append(legend) 
      c_11.Update() 


  c_0 .SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsmu_variations.pdf")
  c_1 .SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dstau_variations.pdf")
  c_10.SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsstarmu_variations.pdf")
  c_11.SaveAs(f"{dest}/{var_bin[0]}_in_bin_{var_bin[1]}_dsstartau_variations.pdf")

 

  
