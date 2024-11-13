import matplotlib as plt
import numpy as np
import ROOT
from ROOT import TLorentzVector
from ROOT import TVector3
import sys
import os
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

from cms_style import CMS_lumi
from copy import deepcopy as dc

# PERFROM SIGNFLIP FIT 
# FIND: B,S


#############################################
#                                           #
#                                           #
#                                           #
#                   *                       #
#                 *   *                     #
#                 *   *                     #
#                *     *                    #
#               *       *                   #
#              *         *                  #
#             *            *                #
#            *       S       *              #
#          *___________________*            #
#       *                          *        #
# *                  B                   *  #
#                                           #
#############################################  Ds mass



############## FIT TO MC ################
 
# Binning:
bins  = 21
start = 1.91#1.91
stop  = 2.028#2.028

#plotting settings
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True) #important for error calc.!
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadLeftMargin(0.15)


outdir = "/work/pahwagne/RDsTools/comb/signflip_fits/"

#this command is a test for git

def getSignflipRatio(hBkg_kk, hBkg_pimu, hSig, hData):

  print(" ======> Prepare variables and fit")
  #path to save
  os.system(f"mkdir -p {outdir}/plots/")


  # initial values
  n0_bkg_kk      = hBkg_kk.Integral()
  n0_bkg_pimu    = hBkg_pimu.Integral()
  n0_sig    = hSig.Integral()
  n0_Data   = hData.Integral() 

 
  #scale signal to comb
  #hSig.Scale(n0_bkg / n0_sig)

  # and create ptr to this histo ( for the cosntructor of hRoo)
  #hData_Ptr = hData.GetPtr()
  #hBkg_Ptr  = hBkg.GetPtr()
  #hSig_Ptr  = hSig.GetPtr()

  var = "phiPi_m"
  # Fitting variable
  dsMass = ROOT.RooRealVar(var,r"D_{s} mass", start, stop)
  dsMass.setRange("complete" ,start ,stop)

  # Define a RooDataHist from TH1D
  hData_Roo      = ROOT.RooDataHist("hData_Roo", "hData_Roo", ROOT.RooArgList(dsMass), hData)
  #hBkg_Roo      = ROOT.RooDataHist("hBkg_Roo" , "hBkg_Roo" , ROOT.RooArgList(dsMass), hBkg)
  hBkg_kk_Roo    = ROOT.RooDataHist("hBkg_kk_Roo" , "hBkg_kk_Roo" , ROOT.RooArgList(dsMass), hBkg_kk)
  hBkg_pimu_Roo  = ROOT.RooDataHist("hBkg_pimu_Roo" , "hBkg_pimu_Roo" , ROOT.RooArgList(dsMass), hBkg_pimu)
  hSig_Roo       = ROOT.RooDataHist("hSig_Roo" , "hSig_Roo" , ROOT.RooArgList(dsMass), hSig)

  #create pdf from histo
  #pdf_bkg     = ROOT.RooHistPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgSet(dsMass), hBkg_Roo);
  pdf_bkg_kk     = ROOT.RooHistPdf("pdf_kk_bkg", "pdf_kk_bkg", ROOT.RooArgSet(dsMass), hBkg_kk_Roo);
  pdf_bkg_pimu   = ROOT.RooHistPdf("pdf_pimu_bkg", "pdf_pimu_bkg", ROOT.RooArgSet(dsMass), hBkg_pimu_Roo);
  pdf_sig     = ROOT.RooHistPdf("pdf_sig", "pdf_sig", ROOT.RooArgSet(dsMass), hSig_Roo);

  kk_scale   = ROOT.RooRealVar( "kk_scale",     "kk_scale"      ,  0.5   ,    0.01, 01.0)
  pimu_scale  = ROOT.RooFormulaVar("pimu_scale", "1.0 - kk_scale", ROOT.RooArgList(kk_scale))

  bkg_scale   = ROOT.RooRealVar( "bkg_scale",     "bkg_scale",  0.8   ,    0.01, 1.0)
  sig_scale  = ROOT.RooFormulaVar("sig_scale", "1.0 - bkg_scale", ROOT.RooArgList(bkg_scale))

  nData       = ROOT.RooRealVar( "nData",          "nData"   , n0_Data,    n0_Data*0.5, n0_Data*1.5)

  #sig_scale   = ROOT.RooRealVar( "sig_scale",     "sig_scale",  1  ,    0.00001, 10)

  # Define functions

  pdf_bkg       = ROOT.RooAddPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgList(pdf_bkg_kk,pdf_bkg_pimu), ROOT.RooArgList(kk_scale, pimu_scale))
  pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg, pdf_sig), ROOT.RooArgList(bkg_scale, sig_scale ))
  #pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg,pdf_sig), ROOT.RooArgList(bkg_scale))
  pdf_total_ext = ROOT.RooExtendPdf("pdf_total_scaled", "pdf_total_scaled", pdf_total, nData, "complete") 
 
  #perform the fit
  #result = pdf_total.fitTo(hData_Roo, ROOT.RooFit.Range("dsMassRange"), ROOT.RooFit.SumW2Error(True))
  result = pdf_total_ext.fitTo(hData_Roo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  
  #get results of parameters
  k_scale    = round(kk_scale.getVal()  ,5)
  p_scale    = round(pimu_scale.getVal()  ,5)
  b_scale    = round(bkg_scale.getVal()  ,5)
  s_scale    = round(sig_scale.getVal()  ,5)
  n_scale    = round(nData.getVal()  ,5)
  #s_scale = round(sig_scale.getVal()  ,5)

  print(f" =======> kk scale is: {k_scale} ")
  print(f" =======> pimu scale is: {p_scale} ")
  print(f" =======> bkg scale is: {b_scale} ")
  print(f" =======> sig scale is: {s_scale} ")
  print(f" =======> norm scale is: {n_scale} ")

  #plotting
  print(f" =======> Start plotting")

  #create frame
  frame = dsMass.frame(ROOT.RooFit.Title(var))
  frame.GetXaxis().SetTitle("D_{s} mass")
  frame.GetXaxis().SetTitleOffset(1.5)
  frame.GetYaxis().SetRangeUser(0, 200000)


  #first plot data and model, they define the normalization range 
  hData_Roo.plotOn(frame,ROOT.RooFit.Name("data"),  ROOT.RooFit.SumW2Error(True)  )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kYellow -7 ),       ROOT.RooFit.Name("pdf_total"), ROOT.RooFit.FillColor(ROOT.kYellow  - 7), ROOT.RooFit.DrawOption("F") )

  #get reduce chi2, 8 is the number of parameters
  chi2  = round(frame.chiSquare("pdf_total","hData_Roo", 2), 5)

  #plot only RooFit method

  c2 = ROOT.TCanvas("Roo","Roo",700,700)
  pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.89, 0.99)
  pad1.SetLeftMargin(0.15)
  pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.89, 0.2)
  pad2.SetLeftMargin(0.15)
  pad3 = ROOT.TPad("pad3", "pad3", 0.82, 0.01, 0.99, 0.2)
  pad3.SetLeftMargin(0.15)
  
  pad1.Draw()
  pad2.Draw()
  pad3.Draw()
  
  #plot subcomponents
  #pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray),      ROOT.RooFit.FillColor(ROOT.kGray),      ROOT.RooFit.Name("pdf_bkg"),          ROOT.RooFit.Components("pdf_bkg"),         ROOT.RooFit.DrawOption("F"),  ROOT.RooFit.FillStyle(1001))
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.DrawOption("F") , ROOT.RooFit.FillStyle(3444)        )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7)        )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), )

  #hist_bkg = pdf_bkg.createHistogram("hist_bkg", dsMass)
  #hist_sig = pdf_sig.createHistogram("hist_sig", dsMass)

  #plot data again, s.t. it is visible on top
  hData_Roo.plotOn(frame,ROOT.RooFit.Name("hData_Roo"))
  
  #get pulls and add them to frame2
  hpull = frame.pullHist("hData_Roo","pdf_total")
  frame2 = dsMass.frame(ROOT.RooFit.Title(""))
  frame2.addPlotable(hpull,"P")
  
  c2.cd()
  pad1.cd()

  #hist_bkg.Draw("HIST")
  #hData.Draw("HIST SAME")
  #hist_sig.Draw("HIST SAME")
  frame.Draw("SAME")
  #hist_bkg.Draw("HIST SAME")


  leg2 = ROOT.TLegend(0.55,0.7,0.8,0.88);
  leg2.SetFillColor(ROOT.kWhite);
  leg2.SetLineColor(ROOT.kWhite);
  leg2.SetTextSize(0.030)
  leg2.AddEntry("hData_Roo","Data","EP")
  leg2.AddEntry("pdf_bkg","  Wrong Sign Comb.","F")
  leg2.AddEntry("pdf_sig","  Signal + Hb","F")
  leg2.AddEntry("pdf_total","  Total","F")
  leg2.Draw("SAME")
  
  text = ROOT.TPaveText(.18,.7,.5,.88,"brNDC");
  text.SetFillColor(0)
  text.SetFillStyle(0)
  text.SetTextAlign(11)
  leg2.SetTextSize(0.035)
  text.AddText(" Relative Scale = " + f"{b_scale} ")
  text.SetBorderSize(0)
  text.Draw("SAME")
  
  CMS_lumi(pad1, 4, 0, cmsText = '    CMS', extraText = '      Private work', lumi_13TeV = '')
  
  pad2.cd()
  frame2.Draw()
  frame2.GetYaxis().SetTitle("#Delta / #sigma")
  frame2.GetXaxis().SetTitle("")
  frame2.GetXaxis().SetTitleSize(3.* frame2.GetXaxis().GetTitleSize())
  frame2.GetXaxis().SetTickLength(0.1)
  frame2.GetYaxis().SetTitleSize(4.* frame2.GetYaxis().GetTitleSize())
  frame2.GetXaxis().SetLabelSize(3.* frame2.GetXaxis().GetLabelSize())
  frame2.GetYaxis().SetLabelSize(3.* frame2.GetYaxis().GetLabelSize())
  frame2.GetYaxis().SetTitleOffset( 0.34 )
  frame2.GetYaxis().CenterTitle()
  
  #divisons on residual y axis
  ratio_bins = 10
  frame2.GetYaxis().SetNdivisions(ratio_bins)
  #borders of residual plot
  ratio_start = -4 
  ratio_stop = 4
  frame2.GetYaxis().SetRangeUser(ratio_start,ratio_stop)
  
  #plot some lines at +- 1 uncertainty of residuals
  l1 = ROOT.TLine(start,1,stop,1)
  l2 = ROOT.TLine(start,-1,stop,-1)
  l1.SetLineStyle(2)
  l2.SetLineStyle(2)
  l1.Draw()
  l2.Draw()
  
  #fill residual histogram
  pad3.cd()
  hgauss = ROOT.TH1D("gauss","gauss",ratio_bins,ratio_start,ratio_stop)
  
  for i in range(bins):
  	#get value of residual and fil it in hgauss
  	event = hpull.Eval(hData.GetBinCenter(i+1))
  	a = hgauss.Fill(event)
  
  #plot residual distribution
  hgauss.GetYaxis().SetTickSize(0)
  hgauss.GetXaxis().SetTickSize(0.1)
  hgauss.GetXaxis().SetNdivisions(ratio_bins)
  hgauss.GetYaxis().SetLabelSize(0)
  hgauss.GetXaxis().SetLabelSize(3.* hgauss.GetXaxis().GetLabelSize())
  hgauss.SetLineColor(ROOT.kBlack)
  hgauss.Draw("hbar")
  
  c2.Modified()
  c2.Update()
  c2.SaveAs(outdir + "signflip_mass.pdf")

  return k_scale, b_scale, n_scale

