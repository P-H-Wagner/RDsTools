import ROOT
import array
import os
import sys
import numpy as np

#take the same pt binning as for fs :)
from get_fs     import bins_fs as bins_pt
from get_fs     import values_fs 

sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

ROOT.gStyle.SetOptTitle(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() #apply weights!


#create RDF 
#input file is a gen lvl sample wout filter
chain = ROOT.TChain("tree")
chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dsmu_with_y}/*.root")
df = ROOT.RDataFrame(chain)

#define reasonable y bins
bins_y = [
0,
0.25,
0.5,
0.75,
1,
1.25,
1.5,
1.75,
2,
#2.5,
#3,
#3.5,
#4,
#10

]

#now we will do a 2D histogram in (pt,y)


pt_edges = array.array("d", bins_pt)
y_edges  = array.array("d", bins_y )

model = ROOT.RDF.TH2DModel("h2", "h2", len(pt_edges)-1, pt_edges, len(y_edges )-1, y_edges)

#get the pythia histo and normalize
pythia = df.Define("abs_gen_bs_y","abs(gen_bs_y)").Histo2D(model, "gen_bs_pt","abs_gen_bs_y")
pythia.Scale(1.0/pythia.Integral())
pythia.GetXaxis().SetTitle("p_{T}")
pythia.GetYaxis().SetTitle("|y|")
pythia.GetZaxis().SetTitle("d#sigma/dp_{T}/dy (normalized)")
pythia.GetZaxis().SetTitleOffset(2.0)
pythia.SetTitle("PYTHIA")

# similarly, import the x-section values for FONLL
from fonll import histo_fonll

#same binning etc
fonll = pythia.Clone()
fonll.SetTitle("FONLL")

for i in range(1,pythia.GetNbinsX()+1):
  for j in range(1,pythia.GetNbinsY()+1):

    fonll.SetBinContent(i,j,histo_fonll[i-1][j-1] * values_fs[i-1])

#normalize
fonll.Scale(1.0/fonll.Integral())


ratio = pythia.Clone()
ratio.SetTitle("FONLL / PYTHIA")
ratio.GetZaxis().SetTitle("a.u.")

fonll_over_pythia = []

for i in range(1,pythia.GetNbinsX()+1):

  dummy = []
  for j in range(1,pythia.GetNbinsY()+1):

    print(f"Comparing: ({i},{j}): ", pythia.GetBinContent(i,j), fonll.GetBinContent(i,j))
    if pythia.GetBinContent(i,j) != 0.0:
      dummy.append(fonll.GetBinContent(i,j) / pythia.GetBinContent(i,j))
      ratio.SetBinContent(i,j,fonll.GetBinContent(i,j) / pythia.GetBinContent(i,j))
    else:
      dummy.append(np.nan)
      ratio.SetBinContent(i,j,np.nan)

  fonll_over_pythia.append(dummy)

#save weights to json
with open(f"fonll_over_pythia.json","w") as f:
  json.dump(fonll_over_pythia, f)
with open(f"bins_pt.json","w") as f:
  json.dump(bins_pt, f)
with open(f"bins_y.json","w") as f:
  json.dump(bins_y, f)

c1 = ROOT.TCanvas('', '', 700, 700)
c1.SetRightMargin(0.20)
pythia.Draw("COLZ")
c1.SaveAs("pythia.pdf")

c2 = ROOT.TCanvas('', '', 700, 700)
c2.SetRightMargin(0.20)
fonll.Draw("COLZ")
c2.SaveAs("fonll.pdf")

c3 = ROOT.TCanvas('', '', 700, 700)
c3.SetRightMargin(0.20)
ratio.Draw("COLZ")
c3.SaveAs("ratio.pdf")
