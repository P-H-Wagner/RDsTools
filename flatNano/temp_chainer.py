import ROOT
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

#create rdf
files = "HOOK_FILES_IN" # this is a list of strings
chain = ROOT.TChain("Events")
for f in files: chain.Add(f)

#output
destination = "HOOK_FILE_OUT"

print("saving to:", destination)
print("====> Create rdf")
df = ROOT.RDataFrame(chain)
print("====> rdf DONE")

print("====> Define branch and create snapshot")
df.Snapshot("tree", destination)
print("====> Snapshot DONE")
