import ROOT
import array
import os
import sys
import numpy as np

#take the same pt binning as for fs :)
from fs_over_fu import bins as bins_pt

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

#create RDF 
#input file is a gen lvl sample wout filter
chain = ROOT.TChain("tree")
chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dsmu_with_y}/*.root")
df = ROOT.RDataFrame(chain)

#define reasonable y bins
bins_y = [
-10,
-4,
-3.5,
-3,
-2.5,
-2,
-1.5,
-1,
-0.5,
0,
0.5,
1,
1.5,
2,
2.5,
3,
3.5,
4,
10

]

#now we will do a 2D histogram in (pt,y)


pt_edges = array.array("d", bins_pt)
y_edges  = array.array("d", bins_y )

model = ROOT.RDF.TH2DModel("h2", "h2", len(pt_edges)-1, pt_edges, len(y_edges )-1, y_edges)

#get the pythia histo and normalize
pythia = df.Histo2D(model, "gen_bs_pt","bs_y")
pythia.Scale(1.0/pythia.Integral())


# similarly, import the x-section values for FONLL

