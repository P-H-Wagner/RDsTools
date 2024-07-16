import ROOT
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *
import subprocess

#create rdf
files = HOOK_FILE_IN # this is a list of strings


#chain = ROOT.TChain("Events")
#for f in files: chain.Add(f)

#print(f"Chaining the files: \n {files}")

#output
destination = "HOOK_FILE_OUT"

print("saving to:", destination)

result = subprocess.run(['hadd', '-f', destination] + files, check = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True )

#print("====> Create rdf")
#df = ROOT.RDataFrame(chain)
#print("====> rdf DONE")

#print("====> Define branch and create snapshot")
#df.Snapshot("tree", destination)
#print("====> Snapshot DONE")
