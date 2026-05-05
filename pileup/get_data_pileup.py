import pdb
import re
import ROOT
import numpy as np
import argparse
import glob
import sys
import ROOT

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

triggers = [
"Mu7_IP4",
"Mu8_IP3",
"Mu8_IP5",
"Mu8p5_IP3p5",
"Mu9_IP5",
"Mu9_IP6",
"Mu12_IP6",
]

#contains histos, name type: pileup_2018D_HLT_Mu7_IP4_run_from_320673_to_322349_GOLDEN.root
histos = glob.glob("/work/manzoni/CMSSW_10_6_25/src/pileup*GOLDEN.root")
#contains trees, name type: 2018D_HLT_Mu7_IP4_run_from_320673_to_322349.root
trees = glob.glob("/work/manzoni/CMSSW_10_6_25/src/2018*.root")

histos = [f for f in histos if any(t in f for t in triggers)]
trees  = [f for f in trees  if any(t in f for t in triggers)]

#debug
#histos  = histos[:2]

#will hold tuple of histo file and #events
files = []

for h in histos:

  #extract run
  pileup_trigger = h.split("_run_from_")[0]
  trigger        = pileup_trigger.split("pileup_")[1]

  run_GOLDEN     = h.split("_run_from_")[1]
  run            = run_GOLDEN.split("_GOLDEN")[0] #str of type 320673_to_322349

  #find the corresponding tree
  found = 0
  for tr in trees:
    if (trigger in tr) and (run in tr):
      if found > 0: print("ALERT, found two matching trees!"); sys.exit()
      else        : 
        print(f"Matching histo - tree pair found \n {h} - {tr}")
        found +=1
      
        nevents = ROOT.RDataFrame("tree",tr).Count().GetValue()
        files.append((h,nevents,f"{trigger}_run_from_{run}"))

  #files = [
  #    ("/work/manzoni/CMSSW_10_6_25/src/pileup_2018D_HLT_Mu7_IP4_run_from_320673_to_322349_GOLDEN.root", 358254225),
  #    ("/work/manzoni/CMSSW_10_6_25/src/pileup_2018D_HLT_Mu7_IP4_run_from_322349_to_324295_GOLDEN.root", 196614520),
  #    ("/work/manzoni/CMSSW_10_6_25/src/pileup_2018D_HLT_Mu7_IP4_run_from_324295_to_325175_GOLDEN.root", 239917306),
  #]


# we do this for all triggers

# Legend
pattern = r"pileup_2018D_(HLT_[^_]+_[^_]+)_run_from_(\d+)_to_(\d+)_GOLDEN"
replacement = r"\1 runs \2 - \3"

#shorts = [re.sub(pattern, replacement, ff[0]) for ff in files]
shorts  = [ff[2] for ff in files]

denominator = np.sum([ff[1] for ff in files])

# Output ROOT file
#outFile = ROOT.TFile("pileup_2018D_HLT_Mu7_IP4_GOLDEN_combined.root", "RECREATE")
outFile = ROOT.TFile("pileup_GOLDEN_combined.root", "RECREATE")

# Create THStack
stack      = ROOT.THStack("pileup_stack"     , ";N_{PU};Fraction")
stack_incl = ROOT.THStack("pileup_stack_incl", ";N_{PU};Fraction")


colors = [
  ROOT.kCyan, 
  ROOT.kCyan + 2, 
  ROOT.kCyan + 4,
  ROOT.kAzure ,
  ROOT.kViolet ,
  ROOT.kGreen, 
  ROOT.kViolet + 4,
  ROOT.kGreen - 5,
  ROOT.kMagenta,
  ROOT.kMagenta - 5,
  ROOT.kMagenta + 3,
  ROOT.kYellow ,
  ROOT.kOrange + 7,
  ROOT.kOrange + 10,
 

  ]  # for stacked histograms

colors_trigger = {

"Mu7_IP4": ROOT.kYellow,
"Mu8_IP3": ROOT.kCyan,
"Mu8_IP5": ROOT.kCyan + 2,
"Mu8p5_IP3p5": ROOT.kCyan -5,
"Mu9_IP5": ROOT.kMagenta,
"Mu9_IP6": ROOT.kMagenta -5,
"Mu12_IP6": ROOT.kGreen,
}

for t in triggers: colors_trigger

normalized_hists      = []  # keep references


for i, (fname, target_yield, label) in enumerate(files):
    print(i)
    f = ROOT.TFile.Open(fname)
    h = f.Get("pileup")
    if not h:
        raise RuntimeError(f"Histogram 'pileup' not found in {fname}")
        
    h = h.Clone(f"pileup_{fname}")  # clone to avoid file ownership issues
    h.SetDirectory(0)  # detach from input file
    
    # Normalize to target yield
    integral = h.Integral()
    if integral > 0:
        scale = float(target_yield) / integral / denominator
        h.Scale(scale)
    else:
        #raise RuntimeError(f"Histogram in {fname} has zero integral!")
        print(f"Histogram in {fname} has zero integral!")
        continue

    # Style
    h.SetFillColor(colors[i % len(colors)])
    #h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(0)
    
    # Add to stack
    stack.Add(h)

    # Add to total
    if i == 0:
         print("first iteration im here!")
        
         total = h.Clone()
         total.SetLineColor(ROOT.kBlack)
         total.SetLineWidth(2)
         total.SetDirectory(0)
    else:
         total.Add(h)

    # Keep normalized version
    normalized_hists.append(h)


#################################
# Plot all single contributions #
#################################

# Plot stack + total
c = ROOT.TCanvas("c", "", 700, 700)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)
stack.Draw("HIST")
#total = stack.GetHistogram()
#for hh in normalized_hists:
#    total.Add(hh)
total.SetFillStyle(0)
total.SetLineWidth(3)
total.SetLineColor(ROOT.kBlack)
total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")

leg = ROOT.TLegend(0.4, 0.4, 0.88, 0.88)
leg.SetBorderSize(0)
for i, h in enumerate(normalized_hists):
    leg.AddEntry(h, f"{shorts[i]} <#mu> = {h.GetMean():.1f}", "f")
leg.AddEntry(total, f"Total <#mu> = {total.GetMean():.1f}", "l")
leg.Draw('SAME')

c.SaveAs("pileup_stack_lin.pdf")

ROOT.gPad.SetLogy(True)

total.Draw("HIST")
stack.Draw("HIST")
total.Draw("HIST SAME")
leg.Draw('SAME')
ROOT.gPad.Modified()
ROOT.gPad.Update()
c.SaveAs("pileup_stack_log.pdf")


print("✅ Histograms, stack, and plot saved successfully!")

#for ibin in range(total.GetNbinsX()):
#    print(f"bin {ibin+1} \t content {total.GetBinContent(ibin+1)}")


# Save everything to output ROOT file
outFile.cd()
for h in normalized_hists:
    h.Write()         # save individual normalized histograms
# stack.Write()         # save stack
total.SetName('total')
total.Write()         # save total
outFile.Close()


############################################################

stack      = ROOT.THStack("pileup_stack"     , ";N_{PU};Fraction")
normalized_hists      = []  # keep references

shorts_reordered = []

for t in triggers:
  
  for i, (fname, target_yield, label) in enumerate(files):
     
      if t not in label: continue

      shorts_reordered.append(label)

      f = ROOT.TFile.Open(fname)
      h = f.Get("pileup")
      if not h:
          raise RuntimeError(f"Histogram 'pileup' not found in {fname}")
          
      h = h.Clone(f"pileup_{fname}")  # clone to avoid file ownership issues
      h.SetDirectory(0)  # detach from input file
      
      # Normalize to target yield
      integral = h.Integral()
      if integral > 0:
          scale = float(target_yield) / integral / denominator
          h.Scale(scale)
      else:
          #raise RuntimeError(f"Histogram in {fname} has zero integral!")
          print(f"Histogram in {fname} has zero integral!")
          continue
  
      # Style
      h.SetFillColor(colors_trigger[t])
      #h.SetLineColor(ROOT.kBlack)
      h.SetLineWidth(0)
      
      # Add to stack
      stack.Add(h)
  
   
      # Add to total
      if i == 0:
           print("first iteration im here!")
          
           total = h.Clone()
           total.SetLineColor(ROOT.kBlack)
           total.SetLineWidth(2)
           total.SetDirectory(0)
      else:
           total.Add(h)
  
      # Keep normalized version
      normalized_hists.append(h)


#################################
# Plot ordered after trigger    # 
#################################

# Plot stack + total
c = ROOT.TCanvas("c", "", 700, 700)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)
stack.Draw("HIST")
total.SetFillStyle(0)
total.SetLineWidth(3)
total.SetLineColor(ROOT.kBlack)
total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")

leg = ROOT.TLegend(0.4, 0.4, 0.88, 0.88)
leg.SetBorderSize(0)
for i, h in enumerate(normalized_hists):
    leg.AddEntry(h, f"{shorts_reordered[3][i]} <#mu> = {h.GetMean():.1f}", "f")
leg.AddEntry(total, f"Total <#mu> = {total.GetMean():.1f}", "l")
leg.Draw('SAME')

c.SaveAs("pileup_stack_incl_lin.pdf")

ROOT.gPad.SetLogy(True)

total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")
leg.Draw('SAME')
ROOT.gPad.Modified()
ROOT.gPad.Update()
c.SaveAs("pileup_stack_incl_log.pdf")


print("✅ Histograms, stack_incl, and plot saved successfully!")

#for ibin in range(total.GetNbinsX()):
#    print(f"bin {ibin+1} \t content {total.GetBinContent(ibin+1)}")



