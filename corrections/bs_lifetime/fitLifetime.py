import ROOT
import numpy as np
import argparse
from datetime import datetime
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

#create RDF 
#input file is a gen lvl sample wout filter
chain = ROOT.TChain("tree")
chain.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dsmu_with_time}/*")
df = ROOT.RDataFrame(chain)

#variable and bins
var = "time"
start = 0.0
stop  = 1e-11
bins  = 50

#create a histo model
model = ROOT.RDF.TH1DModel(var, "", bins, start, stop)

#create a histo
h     = df.Histo1D(model,"time")
h.SetMarkerStyle(20)
h.SetLineColor(ROOT.kBlack)

#prepare Roofit

#variable domain
time = ROOT.RooRealVar(var, "t", start, stop)
time.setRange("complete", start,stop)

h_Roo = ROOT.RooDataHist("h_Roo", "h_Roo", ROOT.RooArgList(time), h.GetPtr())

#define pdf

nominal = -1.5e12
tau_inv = ROOT.RooRealVar("tau_inv", "tau_inv", nominal, -4e12, -0.5e12 )
N       = ROOT.RooRealVar("N"      , "N"      , 50000  , 0      , 100000)

exp     = ROOT.RooExponential("exp"    , "exp"    , time, tau_inv)
exp_ext = ROOT.RooExtendPdf  ("exp_ext", "exp_ext", exp , N      , "complete")

# fit!
result  = exp_ext.fitTo(h_Roo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))


#results
print(f"====> tau = {abs(1/tau_inv.getVal())}")
print(f"====> N   = {N.getVal()}        ")

# plotting
c = ROOT.TCanvas("Roo","Roo",700,700)
c.SetLeftMargin(0.15) #increase to not cut of y axis label

# define a frame
frame = time.frame()
frame.GetXaxis().SetTitle(r"B_{s} decay time [s]")
frame.SetTitle("")

# plot data
h_Roo.plotOn(frame, 
    ROOT.RooFit.MarkerStyle(20)
)

# plot model
exp_ext.plotOn(frame,
    ROOT.RooFit.Range("complete"),
    ROOT.RooFit.NormRange("complete"),
    ROOT.RooFit.LineColor(ROOT.kRed)
)

# draw frame
frame.Draw()

#legend
leg = ROOT.TLegend(0.7,0.8,0.8,0.85);
leg.SetFillColor(ROOT.kWhite);
leg.SetLineColor(ROOT.kWhite);
leg.SetTextSize(0.030)
leg.AddEntry(h.GetPtr(),r"GEN lvl MC","EP")
leg.Draw("SAME")


#text
text = ROOT.TPaveText(.18,.7,.5,.85,"brNDC");
text.SetFillColor(0)
text.SetFillStyle(0)
text.SetTextAlign(11)
text.AddText(rf" #tau = {round( 1e12 * abs(1/tau_inv.getVal()),4)} e-12 s")
#text.AddText(rf" N    = {round(N.getVal(),4)}")
text.SetBorderSize(0)
text.Draw("SAME")


c.SaveAs("fitted_tau.pdf")
