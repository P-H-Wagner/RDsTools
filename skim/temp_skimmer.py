import ROOT
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

bsMass_ = 5.36688

selections = {"baseline": baseline, "bkg": bkg, "ma_cut": ma_cut, "ma_cut_wout_fv":ma_cut_wout_fv}

# pick the skim selection
name = "HOOK_SELECTION"
selec = selections[name]
# skim it!

#create rdf
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/HOOK_DATE_TIME/*"
chain = ROOT.TChain("tree")
chain.Add(files)

#output
destination = f"/scratch/pahwagne/skimmed_HOOK_SELECTION_HOOK_DATE_TIME.root"
print("saving to:", destination)
print("====> Create rdf")
df = ROOT.RDataFrame(chain)
print("====> rdf DONE")

print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
