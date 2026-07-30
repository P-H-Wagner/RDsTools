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
ROOT.TH1.SetDefaultSumw2();


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
histos = glob.glob("/work/manzoni/CMSSW_10_6_25/src/pileup_2018*run_from*GOLDEN.root")
#contains trees, name type: 2018D_HLT_Mu7_IP4_run_from_320673_to_322349.root
trees = glob.glob("/work/manzoni/CMSSW_10_6_25/src/2018*run_from*.root")

histos = [f for f in histos if any(t in f for t in triggers)]
trees  = [f for f in trees  if any(t in f for t in triggers)]

#debug
#histos  = histos[:3]

#will hold tuple of histo file and #events
files = []

for h in histos:

  #extract trigger
  pileup_trigger = h.split("_run_from_")[0]
  year_trigger   = pileup_trigger.split("pileup_")[1]
  trigger        = year_trigger[6:]

  print(f"Found trigger {trigger}!")

  #extract run
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

        print(f"Number of events are: {nevents}")

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

total.GetXaxis().SetTitle("Pileup")
total.GetYaxis().SetTitle("a.u.")

#save original range
start = total.GetXaxis().GetFirst()
stop  = total.GetXaxis().GetLast ()

total.GetXaxis().SetRange(1, 100)

total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")

leg = ROOT.TLegend(0.4, 0.4, 0.88, 0.88)
leg.SetBorderSize(0)
for i, h in enumerate(normalized_hists):
    leg.AddEntry(h, f"{shorts[i]} <#mu> = {h.GetMean():.1f}", "f")
leg.AddEntry(total, f"Total <#mu> = {total.GetMean():.1f}", "l")
leg.Draw('SAME')

c.SaveAs("pileup_stack_single_runs.pdf")

ROOT.gPad.SetLogy(True)

total.Draw("HIST")
stack.Draw("HIST")
total.Draw("HIST SAME")
leg.Draw('SAME')
ROOT.gPad.Modified()
ROOT.gPad.Update()
c.SaveAs("pileup_stack_single_runs_log.pdf")


print("✅ Histograms, stack, and plot saved successfully!")

#change again on full range for saving!
total.GetXaxis().SetRange(start,stop)

#FILL UP DOWN VARIATION
total_up   = total.Clone("total_up")
total_down = total.Clone("total_down")

for i in range(1,total.GetNbinsX() + 1):

  total_up  .SetBinContent(i, total.GetBinContent(i) + total.GetBinError(i) )
  total_down.SetBinContent(i, total.GetBinContent(i) - total.GetBinError(i) )


# Save everything to output ROOT file
outFile = ROOT.TFile("pileup_GOLDEN_combined.root", "RECREATE")
outFile.cd()

for i,h in enumerate(normalized_hists):
    h.SetName(shorts[i])
    h.Write()         # save individual normalized histograms

# stack.Write()         # save stack
total.SetName('total')

total_up  .Write()
total_down.Write()
total.Write()         # save total

outFile.Close()


############################################################

#############################################
# Plot all contributions grouped by trigger #
#############################################

stack = ROOT.THStack("pileup_stack", ";N_{PU};Fraction")

histos_by_trigger = {}
total = None

for t in triggers:

    this_trigger = 0

    for fname, target_yield, label in files:

        if t not in label:
            continue
       

        print(f"Adding {label} to trigger {t}")

        f = ROOT.TFile.Open(fname)
        h = f.Get("pileup")

        if not h:
            raise RuntimeError(f"Histogram 'pileup' not found in {fname}")

        h = h.Clone()
        h.SetDirectory(0)

        integral = h.Integral()

        if integral <= 0:
            print(f"Histogram in {fname} has zero integral!")
            continue

        scale = float(target_yield) / integral / denominator
        h.Scale(scale)

        if this_trigger == 0:
            combined = h.Clone(f"{t}_combined")
            combined.SetDirectory(0)
        else:
            combined.Add(h)

        #we found a histogram belonging to this trigger
        this_trigger+=1
 
 
    print(f"We found {this_trigger} histograms for this trigger.")
    if this_trigger == 0:
        continue

    combined.SetFillColor(colors_trigger[t])
    combined.SetLineWidth(0)

    histos_by_trigger[t] = combined

    stack.Add(combined)

    if total is None:
        total = combined.Clone("total")
        total.SetDirectory(0)
        total.SetFillStyle(0)
        total.SetLineColor(ROOT.kBlack)
        total.SetLineWidth(3)
    else:
        total.Add(combined)

leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
leg.SetBorderSize(0)

for t in triggers:

    if t not in histos_by_trigger:
        continue

    h = histos_by_trigger[t]

    leg.AddEntry(
        h,
        f"{t}, <#mu> = {h.GetMean():.1f}",
        "f"
    )

leg.AddEntry(total, f"Total <#mu> = {total.GetMean():.1f}", "l")

#################################
# Plot ordered after trigger    #
#################################

c = ROOT.TCanvas("c_incl", "", 700, 700)
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)

stack.Draw("HIST")

total.SetFillStyle(0)
total.SetLineWidth(3)
total.SetLineColor(ROOT.kBlack)

#FILL UP DOWN VARIATION
total_up   = total.Clone("total_up")
total_up.SetLineColor(ROOT.kRed)

total_down = total.Clone("total_down")
total_down.SetLineColor(ROOT.kBlue)

for i in range(1,total.GetNbinsX() + 1):

  total_up  .SetBinContent(i, total.GetBinContent(i) + total.GetBinError(i) )
  total_down.SetBinContent(i, total.GetBinContent(i) - total.GetBinError(i) )

total.GetXaxis().SetTitle("Pileup")
total.GetYaxis().SetTitle("a.u.")

#save original range
start = total.GetXaxis().GetFirst()
stop  = total.GetXaxis().GetLast ()

total     .GetXaxis().SetRange(1, 100)
total_up  .GetXaxis().SetRange(1, 100)
total_down.GetXaxis().SetRange(1, 100)

total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")

leg.Draw("SAME")

c.SaveAs("pileup_stack_by_trigger.pdf")

ROOT.gPad.SetLogy(True)

total.Draw("HIST")
stack.Draw("HIST SAME")
total.Draw("HIST SAME")
leg.Draw("SAME")

ROOT.gPad.Modified()
ROOT.gPad.Update()

c.SaveAs("pileup_stack_by_trigger_log.pdf")


ROOT.gPad.SetLogy(False)

total     .SetLineWidth(1)
total_up  .SetLineWidth(1)
total_down.SetLineWidth(1)
leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
leg.AddEntry(total     , "Pileup central", "L")
leg.AddEntry(total_up  , "Pileup up"     , "L")
leg.AddEntry(total_down, "Pileup down"   , "L")

total     .Draw("HIST")
total_up  .Draw("HIST SAME")
total_down.Draw("HIST SAME")
leg.Draw("SAME")

ROOT.gPad.Modified()
ROOT.gPad.Update()

c.SaveAs("pileup_up_down.pdf")

print("✅ Histograms, stack, and plot saved successfully!")

#change again on full range for saving!
total     .GetXaxis().SetRange(start,stop)
total_up  .GetXaxis().SetRange(start,stop)
total_down.GetXaxis().SetRange(start,stop)


#################################
# Save histograms               #
#################################

outFile2 = ROOT.TFile("pileup_GOLDEN_combined_by_trigger.root","RECREATE")

outFile2.cd()
for t in triggers:

    if t not in histos_by_trigger:
        continue

    h = histos_by_trigger[t]

    h.SetName(f"pileup_{t}")
    h.Write()

total.SetName("total")
total.Write()
total_up  .Write()
total_down.Write()

stack.Write("stack")

outFile2.Close()

print("✅ ROOT file saved successfully!")
