import ROOT
import os
import sys


selec= "HOOK_SELECTION"

# skim it!

#create rdf
files = "HOOK_FILE_IN"
chain = ROOT.TChain("tree")
chain.Add(files)

#output
destination = "HOOK_FILE_OUT"
print("saving to:", destination)
print("====> Create rdf")
df = ROOT.RDataFrame(chain)
print("====> rdf DONE")

print("====> Define branch and create snapshot")
"HOOK_NEW_BRANCH"

print("====> Snapshot DONE")
