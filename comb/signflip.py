import matplotlib as plt
import numpy as np
import ROOT
from ROOT import TLorentzVector
from ROOT import TVector3
import sys
import os
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *
from ROOT import RooMsgService, RooFit
from cms_style import CMS_lumi
from copy import deepcopy as dc
from datetime import datetime

RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
RooMsgService.instance().setStreamStatus(0, False)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;") 

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

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S") 

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



#this command is a test for git

def getSignflipRatio(hBkg_kk, hBkg_pimu,hSig, hData, mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, findcut = None, cut = None):

  if findcut:
    outdir = f"/work/pahwagne/RDsTools/comb/signflip_fits/cut_fits/{cut}/"
  else: 
    outdir = f"/work/pahwagne/RDsTools/comb/signflip_fits/{dt}/"

  if not os.path.exists(outdir):
    os.makedirs(outdir)

  print(" ======> Prepare variables and fit")


  # initial values
  n0_bkg_kk      = hBkg_kk.Integral()
  n0_bkg_pimu    = hBkg_pimu.Integral()
  #n0_bkg_both    = hBkg_both.Integral()
  n0_sig    = hSig.Integral()
  n0_Data   = hData.Integral() 

  print("prefit nunmber of events inside function of kk wrong", n0_bkg_pimu)
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
  #hBkg_both_Roo  = ROOT.RooDataHist("hBkg_both_Roo" , "hBkg_both_Roo" , ROOT.RooArgList(dsMass), hBkg_both)
  hSig_Roo       = ROOT.RooDataHist("hSig_Roo" , "hSig_Roo" , ROOT.RooArgList(dsMass), hSig)

  #create pdf from histo
  #pdf_bkg     = ROOT.RooHistPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgSet(dsMass), hBkg_Roo);
  pdf_bkg_kk     = ROOT.RooHistPdf("pdf_kk_bkg", "pdf_kk_bkg", ROOT.RooArgSet(dsMass), hBkg_kk_Roo);
  pdf_bkg_pimu   = ROOT.RooHistPdf("pdf_pimu_bkg", "pdf_pimu_bkg", ROOT.RooArgSet(dsMass), hBkg_pimu_Roo);
  #pdf_bkg_both   = ROOT.RooHistPdf("pdf_both_bkg", "pdf_both_bkg", ROOT.RooArgSet(dsMass), hBkg_both_Roo);
  pdf_sig     = ROOT.RooHistPdf("pdf_sig", "pdf_sig", ROOT.RooArgSet(dsMass), hSig_Roo);

  kk_scale   = ROOT.RooRealVar( "kk_scale",   "kk_scale"      ,  0.5   ,    0.00, 1.0)
  pimu_scale  = ROOT.RooFormulaVar("pimu_scale", "1.0 - kk_scale"    ,  ROOT.RooArgList(kk_scale))

  kk_scale     = ROOT.RooRealVar( "kk_scale",   "kk_scale"      ,  0.0 )
  pimu_scale   = ROOT.RooRealVar( "pimu_scale",   "pimu_scale"      ,  1.0 )


  #pimu_scale  = ROOT.RooRealVar("pimu_scale", "pimu_scale"    ,  0.3   ,    0.00, 1.0)
  #both_scale  = ROOT.RooFormulaVar("both_scale", "1.0 - kk_scale - pimu_scale", ROOT.RooArgList(kk_scale, pimu_scale))

  bkg_scale   = ROOT.RooRealVar( "bkg_scale",     "bkg_scale",  0.5  ,    0.0, 1.0)
  sig_scale  = ROOT.RooFormulaVar("sig_scale", "1.0 - bkg_scale", ROOT.RooArgList(bkg_scale))

  nData       = ROOT.RooRealVar( "nData",          "nData"   , n0_Data,    n0_Data*0.5, n0_Data*1.5)

  #sig_scale   = ROOT.RooRealVar( "sig_scale",     "sig_scale",  1  ,    0.00001, 10)

  # Define functions

  pdf_bkg       = ROOT.RooAddPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgList(pdf_bkg_kk,pdf_bkg_pimu), ROOT.RooArgList(kk_scale, pimu_scale))
  #pdf_bkg       = ROOT.RooAddPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgList(pdf_bkg_kk,pdf_bkg_pimu, pdf_bkg_both), ROOT.RooArgList(kk_scale, pimu_scale, both_scale))

  if n0_sig > 0:
    pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg, pdf_sig), ROOT.RooArgList(bkg_scale, sig_scale ))
  else:
    pdf_total     = pdf_bkg #ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg), ROOT.RooArgList(bkg_scale))

  #pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg,pdf_sig), ROOT.RooArgList(bkg_scale))
  pdf_total_ext = ROOT.RooExtendPdf("pdf_total_scaled", "pdf_total_scaled", pdf_total, nData, "complete") 
 
  #perform the fit
  #result = pdf_total.fitTo(hData_Roo, ROOT.RooFit.Range("dsMassRange"), ROOT.RooFit.SumW2Error(True))
  result = pdf_total_ext.fitTo(hData_Roo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  
  #get results of parameters
  k_scale    = round(kk_scale.getVal()  ,5)
  p_scale    = round(pimu_scale.getVal()  ,5)
  #kp_scale   = round(both_scale.getVal()  ,5)
  b_scale    = round(bkg_scale.getVal()  ,5)
  s_scale    = round(sig_scale.getVal()  ,5)
  n_data     = round(nData.getVal()  ,5)
  n_sig      = round(nData.getVal() * (1 - bkg_scale.getVal())  ,5)
  n_bkg      = round(nData.getVal() * (bkg_scale.getVal())      ,5)
  n_kk       = round(nData.getVal() * (bkg_scale.getVal()) * kk_scale.getVal()       ,5)
  n_pimu     = round(nData.getVal() * (bkg_scale.getVal()) * pimu_scale.getVal()     ,5)
  #s_scale = round(sig_scale.getVal()  ,5)


  print(f" =======> kk scale is: {k_scale} ")
  print(f" =======> pimu scale is: {p_scale} ")
  #print(f" =======> both scale is: {kp_scale} ")
  print(f" =======> bkg scale is: {b_scale} ")
  print(f" =======> sig scale is: {s_scale} ")
  print(f" =======> number of events are: {n_data} ")
  print(f" =======> number of bkg events: {n_data * b_scale}")
  print(f" =======> number of kk bkg events: {n_pimu}")
  print(f" =======> number of pimu bkg events: {n_kk}")
  print(f" =======> number of signal events: {n_data * (1-b_scale)}")

  #plotting
  print(f" =======> Start plotting")

  #create frame
  frame = dsMass.frame(ROOT.RooFit.Title(var))
  frame.GetXaxis().SetTitle("D_{s} mass")
  frame.GetXaxis().SetTitleOffset(1.5)
  #frame.GetYaxis().SetRangeUser(0, 200000)
  max_y = hData.GetMaximum() * 1.2  # Scale the max value slightly for better visualization
  frame.GetYaxis().SetRangeUser(0, max_y)


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
  #pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen -7),   ROOT.RooFit.Name("pdf_kk_bkg"),          ROOT.RooFit.Components("pdf_kk_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen -7), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  #pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen -7),   ROOT.RooFit.Name("pdf_kk_bkg"),          ROOT.RooFit.Components("pdf_kk_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen -7), )
  #pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen +9),   ROOT.RooFit.Name("pdf_pimu_bkg"),          ROOT.RooFit.Components("pdf_pimu_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen +9), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  #pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen +9),   ROOT.RooFit.Name("pdf_pimu_bkg"),          ROOT.RooFit.Components("pdf_pimu_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen +9), )

  hist_bkg = pdf_bkg.createHistogram("hist_bkg", dsMass)
  hist_sig = pdf_sig.createHistogram("hist_sig", dsMass)


  
  print(f"After the fit we have: {hist_sig.Integral()} events")


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
  leg2.AddEntry("pdf_bkg","  Wrong Sign (Comb.)","F")
  #leg2.AddEntry("pdf_kk_bkg"," kk wrong","F")
  #leg2.AddEntry("pdf_pimu_bkg"," pimu wrong","F")
  leg2.AddEntry("pdf_sig","  Signal + Hb","F")
  leg2.AddEntry("pdf_total",r" PDF_{ext}","F")
  leg2.Draw("SAME")
  
  text = ROOT.TPaveText(.18,.7,.5,.88,"brNDC");
  text.SetFillColor(0)
  text.SetFillStyle(0)
  text.SetTextAlign(11)
  leg2.SetTextSize(0.035)
  text.AddText(" Relative Scale = " + f"{b_scale} ")
  #text.AddText(" kk Scale = "       + f"{k_scale} ")
  #text.AddText(" pi mu Scale = "    + f"{p_scale} ")
  #text.AddText(" both Scale = "     + f"{kp_scale} ")
  text.SetBorderSize(0)
  text.Draw("SAME")
  
  #CMS_lumi(pad1, 4, 0, cmsText = '    CMS', extraText = '      Private work', lumi_13TeV = '')
  
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

  hist_bkg.Scale(n_bkg/ hist_bkg.Integral())
  #set bin errors manually, what is root doing here!?
  for i in range(1, hist_bkg.GetNbinsX() + 1):  # Bins start at 1 in ROOT
    bin_content = hist_bkg.GetBinContent(i)
    hist_bkg.SetBinError(i, np.sqrt(bin_content) )

  #####################################
  # now fit the sf bkg to a chebyshev #
  #####################################
  hBkgRoo   =  ROOT.RooDataHist("hBkg","hBkg",ROOT.RooArgList(dsMass),hist_bkg)

  a1 = ROOT.RooRealVar ("c1", "Coefficient 1", -0.3, -1.0, 1.0);
  a2 = ROOT.RooRealVar ("c2", "Coefficient 2", 0.2, -1.0, 1.0);
  a3 = ROOT.RooRealVar ("c3", "Coefficient 3", -0.1, -1.0, 1.0);

  mean1   = ROOT.RooRealVar(    "mean1",     "mean1", 1.968,  1.92, 2.0)
  sigma1  = ROOT.RooRealVar(    "sigma1",    "sigma1",  0.012,     0.001,  0.050)
  gauss1      = ROOT.RooGaussian(   "gauss1", "gauss1", dsMass,  mean1, sigma1) 

  mean2   = ROOT.RooRealVar(    "mean2",     "mean2", 1.955,  1.92 , 2.0)
  sigma2  = ROOT.RooRealVar(    "sigma2",    "sigma2",  0.008,     0.001,  0.050)
  gauss2  = ROOT.RooGaussian(   "gauss2", "gauss2", dsMass,  mean2, sigma2) 

  mean3   = ROOT.RooRealVar("mean3", "Mean (log-space)",   1.968, 1.92, 2.0)
  sigma3  = ROOT.RooRealVar("sigma3", "Width (log-space)", 0.01, 0.00, 0.05)

  lognormal = ROOT.RooLognormal("lognormalPdf", "Lognormal PDF", dsMass, mean3, sigma3)


  n_bkg_evt  = ROOT.RooRealVar(    "n_bkg_event",    "n_bkg_event", n_bkg)

  fracs = ROOT.RooRealVar(    "fracs",    "fracs",  0.5,     0.0,  1.0)
  fracs2 = ROOT.RooRealVar(    "fracs2",    "fracs2",  0.5,     0.0,  1.0)

  dg = ROOT.RooAddPdf( "dg","dg PDF",ROOT.RooArgList(gauss1,gauss2),fracs)
  gaussnormal = ROOT.RooAddPdf( "gaussnormal","gaussnormalPDF",ROOT.RooArgList(dg,lognormal),fracs2)
  #doublegauss = ROOT.RooExtendPdf("doublegauss","doublegauss PDF", dg,n_bkg_evt, "complete")

  #cheby = ROOT.RooChebychev ("cheby", "cheby", dsMass, ROOT.RooArgList(a1, a2, a3));
  #result2 = cheby.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #result2 = doublegauss.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #result2 = gauss1.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  result2 = gaussnormal.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #cheby.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("cheby"),             ROOT.RooFit.Components("cheby"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #doublegauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("doublegauss"),             ROOT.RooFit.Components("doublegauss"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )

  dsMass.setRange("left"     ,mlow3  ,mlow2)   
  dsMass.setRange("right"    ,mhigh2 ,mhigh3)        
  dsMass.setRange("signal"   ,mlow  ,mhigh)

  A = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"left"  ).getVal()
  C = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"right" ).getVal()
  B = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"signal").getVal()
  ##import pdb
  print(f"integral values: A = {A}, B = {B}, C = {C}")

  #gauss1.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("gauss1"),             ROOT.RooFit.Components("gauss1"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #gaussnormal.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("gaussnormal"),             ROOT.RooFit.Components("gaussnormal"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #create extended bkg hist
  #hist_bkg_ext = hist_bkg.Clone()
  #hist_bkg_ext.Scale(n_bkg)
  #dg_events = doublegauss.getVal()
  #print("dg events", dg_events)

  #create frame3
  frame3 = dsMass.frame(ROOT.RooFit.Title(var))
  frame3.GetXaxis().SetTitle("D_{s} mass")
  frame3.GetXaxis().SetTitleOffset(1.5)
  #frame3.GetYaxis().SetRangeUser(0, 200000)
  #max_y = hBkgR.GetMaximum() * 1.2  # Scale the max value slightly for better visualization
  frame3.GetYaxis().SetRangeUser(0, max_y)


  #first plot data and model, they define the normalization range 
  hBkgRoo.plotOn(frame3,ROOT.RooFit.Name("data2"),  ROOT.RooFit.SumW2Error(True)  )
  gaussnormal.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed -7 ),       ROOT.RooFit.Name("gaussnormal"), ROOT.RooFit.FillColor(ROOT.kRed  - 7), ROOT.RooFit.DrawOption("F") )
  hBkgRoo.plotOn(frame3,ROOT.RooFit.Name("data3"))
  ##get reduce chi2, 8 is the number of parameters
  #chi2  = round(frame3.chiSquare("pdf_total","hData_Roo", 2), 5)

  ##plot only RooFit method

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
  #
  ##plot subcomponents
  ##pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray),      ROOT.RooFit.FillColor(ROOT.kGray),      ROOT.RooFit.Name("pdf_bkg"),          ROOT.RooFit.Components("pdf_bkg"),         ROOT.RooFit.DrawOption("F"),  ROOT.RooFit.FillStyle(1001))
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.DrawOption("F") , ROOT.RooFit.FillStyle(3444)        )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7)        )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), )


  c2.cd()
  pad1.cd()
  frame3.Draw("SAME")

  c2.Modified()
  c2.Update()
  c2.SaveAs(outdir + "signflip_mass_bebi.pdf")




  print("==============> Saved!")
  return n_kk, n_pimu, n_sig, [A,B,C]


def fitAnotherVar( hBkg_both, hSig, hData, mlow, mhigh, mlow2, mhigh2, mlow3, mhigh3, key, start, stop):

  print(" ======> Prepare variables and fit")

  #path to save
  outdir = f"/work/pahwagne/RDsTools/comb/signflip_fits/{dt}/"
  os.system(f"mkdir -p {outdir}/plots/")


  # initial values
  #n0_bkg_kk      = hBkg_kk.Integral()
  #n0_bkg_pimu    = hBkg_pimu.Integral()
  n0_bkg_both    = hBkg_both.Integral()
  n0_sig    = hSig.Integral()
  n0_Data   = hData.Integral() 

  #print("prefit nunmber of events inside function of kk wrong", n0_bkg_pimu)
  #scale signal to comb
  #hSig.Scale(n0_bkg / n0_sig)

  # and create ptr to this histo ( for the cosntructor of hRoo)
  #hData_Ptr = hData.GetPtr()
  #hBkg_Ptr  = hBkg.GetPtr()
  #hSig_Ptr  = hSig.GetPtr()

  var = "anotherVar"
  # Fitting variable
  dsMass = ROOT.RooRealVar(key,key, start, stop)
  dsMass.setRange("complete" ,start ,stop)

  # Define a RooDataHist from TH1D
  hData_Roo      = ROOT.RooDataHist("hData_Roo"     , "hData_Roo"    , ROOT.RooArgList(dsMass), hData     )
  #hBkg_Roo      = ROOT.RooDataHist("hBkg_Roo"      , "hBkg_Roo"     , ROOT.RooArgList(dsMass), hBkg)
  #hBkg_kk_Roo    = ROOT.RooDataHist("hBkg_kk_Roo"   , "hBkg_kk_Roo"  , ROOT.RooArgList(dsMass), hBkg_kk   )
  #hBkg_pimu_Roo  = ROOT.RooDataHist("hBkg_pimu_Roo" , "hBkg_pimu_Roo", ROOT.RooArgList(dsMass), hBkg_pimu )
  hBkg_both_Roo  = ROOT.RooDataHist("hBkg_both_Roo" , "hBkg_both_Roo" , ROOT.RooArgList(dsMass), hBkg_both)
  hSig_Roo       = ROOT.RooDataHist("hSig_Roo"      , "hSig_Roo"     , ROOT.RooArgList(dsMass), hSig      )

  #create pdf from histo
  #pdf_bkg     = ROOT.RooHistPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgSet(dsMass), hBkg_Roo);
  #pdf_bkg_kk     = ROOT.RooHistPdf("pdf_kk_bkg", "pdf_kk_bkg", ROOT.RooArgSet(dsMass), hBkg_kk_Roo);
  #pdf_bkg_pimu   = ROOT.RooHistPdf("pdf_pimu_bkg", "pdf_pimu_bkg", ROOT.RooArgSet(dsMass), hBkg_pimu_Roo);
  pdf_bkg_both   = ROOT.RooHistPdf("pdf_both_bkg", "pdf_both_bkg", ROOT.RooArgSet(dsMass), hBkg_both_Roo);
  pdf_sig        = ROOT.RooHistPdf("pdf_sig", "pdf_sig", ROOT.RooArgSet(dsMass), hSig_Roo);

  #kk_scale   = ROOT.RooRealVar( "kk_scale",   "kk_scale"      ,  0.5  ,    0.0, 1.0)
  #pimu_scale  = ROOT.RooFormulaVar("pimu_scale", "1.0 - kk_scale"    ,  ROOT.RooArgList(kk_scale))

  kk_scale     = ROOT.RooRealVar( "kk_scale",   "kk_scale"      ,  0.0 )
  pimu_scale   = ROOT.RooRealVar( "pimu_scale",   "pimu_scale"      ,  0.0 )
  both_scale   = ROOT.RooRealVar( "pimu_scale",   "pimu_scale"      ,  1.0 )

  #kk_scale.setError(0.001)

  #include also both wrong
  # Two independent free parameters in [0,1]
  #kk_frac   = ROOT.RooRealVar("kk_frac",   "kk_frac",   0.33, 0.0, 1.0)
  #pimu_frac = ROOT.RooRealVar("pimu_frac", "pimu_frac", 0.5, 0.0, 1.0)
  #
  ## The three fractions that are guaranteed >= 0 and sum to 1
  #kk_scale   = ROOT.RooFormulaVar("kk_scale",   "@0", ROOT.RooArgList(kk_frac))
  #pimu_scale = ROOT.RooFormulaVar("pimu_scale", "(1-@0)*@1", ROOT.RooArgList(kk_frac, pimu_frac))
  #both_scale = ROOT.RooFormulaVar("both_scale", "(1-@0)*(1-@1)", ROOT.RooArgList(kk_frac, pimu_frac))


  bkg_scale   = ROOT.RooRealVar   ("bkg_scale",     "bkg_scale",  0.3  ,    0.0, 1.0)
  sig_scale   = ROOT.RooFormulaVar("sig_scale",     "1.0 - bkg_scale", ROOT.RooArgList(bkg_scale))

  nData       = ROOT.RooRealVar( "nData",          "nData"   , n0_Data,    n0_Data*0.5, n0_Data*1.5)


  # Define functions

  #pdf_bkg       = ROOT.RooAddPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgList(pdf_bkg_kk,pdf_bkg_pimu), ROOT.RooArgList(kk_scale, pimu_scale))
  #pdf_bkg       = ROOT.RooAddPdf("pdf_bkg", "pdf_bkg", ROOT.RooArgList(pdf_bkg_kk,pdf_bkg_pimu, pdf_bkg_both), ROOT.RooArgList(kk_scale, pimu_scale, both_scale))
  pdf_bkg       = pdf_bkg_both 

  if n0_sig > 0:
    pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg, pdf_sig), ROOT.RooArgList(bkg_scale, sig_scale ))
  else:
    pdf_total     = pdf_bkg #ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg), ROOT.RooArgList(bkg_scale))

  #pdf_total     = ROOT.RooAddPdf("pdf_total", "pdf_total", ROOT.RooArgList(pdf_bkg,pdf_sig), ROOT.RooArgList(bkg_scale))
  pdf_total_ext = ROOT.RooExtendPdf("pdf_total_scaled", "pdf_total_scaled", pdf_total, nData, "complete") 
 
  #perform the fit
  #result = pdf_total.fitTo(hData_Roo, ROOT.RooFit.Range("dsMassRange"), ROOT.RooFit.SumW2Error(True))
  result = pdf_total_ext.fitTo(hData_Roo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Strategy(2) )
 

  nll = pdf_total_ext.createNLL(hData_Roo, ROOT.RooFit.SumW2Error(True))
  minim = ROOT.RooMinimizer(nll)
  minim.setStrategy(2)
  minim.minimize("Minuit2", "migrad")
  minim.hesse()
  minim.setEps(1e-6)      
  result = minim.save()
  result.Print("v")


 
  #get results of parameters
  k_scale    = round(kk_scale.getVal()  ,5)
  p_scale    = round(pimu_scale.getVal()  ,5)
  kp_scale   = round(both_scale.getVal()  ,5)
  b_scale    = round(bkg_scale.getVal()  ,5)
  s_scale    = round(sig_scale.getVal()  ,5)
  n_data     = round(nData.getVal()  ,5)
  n_sig      = round(nData.getVal() * (1 - bkg_scale.getVal())  ,5)
  n_bkg      = round(nData.getVal() * (bkg_scale.getVal())      ,5)
  n_kk       = round(nData.getVal() * (bkg_scale.getVal()) * kk_scale.getVal()       ,5)
  n_pimu     = round(nData.getVal() * (bkg_scale.getVal()) * pimu_scale.getVal()     ,5)
  #n_both     = round(nData.getVal() * (bkg_scale.getVal()) * both_scale.getVal()     ,5)
  #s_scale = round(sig_scale.getVal()  ,5)


  print(f" =======> kk scale is: {k_scale} ")
  print(f" =======> pimu scale is: {p_scale} ")
  print(f" =======> both scale is: {kp_scale} ")
  print(f" =======> bkg scale is: {b_scale} ")
  print(f" =======> sig scale is: {s_scale} ")
  print(f" =======> number of events are: {n_data} ")
  print(f" =======> number of bkg events: {n_data * b_scale}")
  print(f" =======> number of kk bkg events: {n_pimu}")
  print(f" =======> number of pimu bkg events: {n_kk}")
  print(f" =======> number of signal events: {n_data * (1-b_scale)}")

  #plotting
  print(f" =======> Start plotting")

  #create frame
  frame = dsMass.frame(ROOT.RooFit.Title(var))
  frame.GetXaxis().SetTitle(key)
  frame.GetXaxis().SetTitleOffset(1.5)
  #frame.GetYaxis().SetRangeUser(0, 200000)
  max_y = hData.GetMaximum() * 1.2  # Scale the max value slightly for better visualization
  frame.GetYaxis().SetRangeUser(0, max_y)


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
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen -7),   ROOT.RooFit.Name("pdf_kk_bkg"),          ROOT.RooFit.Components("pdf_kk_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen -7), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen -7),   ROOT.RooFit.Name("pdf_kk_bkg"),          ROOT.RooFit.Components("pdf_kk_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen -7), )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen +9),   ROOT.RooFit.Name("pdf_pimu_bkg"),          ROOT.RooFit.Components("pdf_pimu_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen +9), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen +9),   ROOT.RooFit.Name("pdf_pimu_bkg"),          ROOT.RooFit.Components("pdf_pimu_bkg"),  ROOT.RooFit.FillColor(ROOT.kGreen +9), )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue +2),   ROOT.RooFit.Name("pdf_both_bkg"),          ROOT.RooFit.Components("pdf_both_bkg"),  ROOT.RooFit.FillColor(ROOT.kBlue + 2), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)  )
  pdf_total.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue +2),   ROOT.RooFit.Name("pdf_both_bkg"),          ROOT.RooFit.Components("pdf_both_bkg"),  ROOT.RooFit.FillColor(ROOT.kBlue + 2), )

  hist_bkg = pdf_bkg.createHistogram("hist_bkg", dsMass)
  hist_sig = pdf_sig.createHistogram("hist_sig", dsMass)


  
  print(f"After the fit we have: {hist_sig.Integral()} events")


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
  leg2.AddEntry("pdf_kk_bkg"," kk wrong","F")
  leg2.AddEntry("pdf_pimu_bkg"," pimu wrong","F")
  #leg2.AddEntry("pdf_both_bkg"," both wrong","F")
  leg2.AddEntry("pdf_sig","  Signal + Hb","F")
  leg2.AddEntry("pdf_total","  Total","F")
  leg2.Draw("SAME")
  
  text = ROOT.TPaveText(.18,.7,.5,.88,"brNDC");
  text.SetFillColor(0)
  text.SetFillStyle(0)
  text.SetTextAlign(11)
  leg2.SetTextSize(0.035)
  text.AddText(" Relative Scale = " + f"{b_scale} ")
  text.AddText(" kk Scale = "       + f"{k_scale} ")
  text.AddText(" pi mu Scale = "    + f"{p_scale} ")
  #text.AddText(" both Scale = "     + f"{kp_scale} ")
  #text.AddText(" both Scale = "     + f"{kp_scale} ")
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
  c2.SaveAs(outdir + f"signflip_{key}.pdf")

  #return n_bkg, n_sig

  hist_bkg.Scale(n_bkg/ hist_bkg.Integral())
  #set bin errors manually, what is root doing here!?
  for i in range(1, hist_bkg.GetNbinsX() + 1):  # Bins start at 1 in ROOT
    bin_content = hist_bkg.GetBinContent(i)
    hist_bkg.SetBinError(i, np.sqrt(bin_content) )

  ######################################
  ## now fit the sf bkg to a chebyshev #
  ######################################
  hBkgRoo   =  ROOT.RooDataHist("hBkg","hBkg",ROOT.RooArgList(dsMass),hist_bkg)

  a1 = ROOT.RooRealVar ("c1", "Coefficient 1", -0.3, -1.0, 1.0);
  a2 = ROOT.RooRealVar ("c2", "Coefficient 2", 0.2, -1.0, 1.0);
  a3 = ROOT.RooRealVar ("c3", "Coefficient 3", -0.1, -1.0, 1.0);

  mean1   = ROOT.RooRealVar(    "mean1",     "mean1", 1.968,  1.92, 2.0)
  sigma1  = ROOT.RooRealVar(    "sigma1",    "sigma1",  0.012,     0.001,  0.050)
  gauss1      = ROOT.RooGaussian(   "gauss1", "gauss1", dsMass,  mean1, sigma1) 

  mean2   = ROOT.RooRealVar(    "mean2",     "mean2", 1.955,  1.92 , 2.0)
  sigma2  = ROOT.RooRealVar(    "sigma2",    "sigma2",  0.008,     0.001,  0.050)
  gauss2  = ROOT.RooGaussian(   "gauss2", "gauss2", dsMass,  mean2, sigma2) 

  mean3   = ROOT.RooRealVar("mean3", "Mean (log-space)",   1.968, 1.92, 2.0)
  sigma3  = ROOT.RooRealVar("sigma3", "Width (log-space)", 0.01, 0.00, 0.05)

  lognormal = ROOT.RooLognormal("lognormalPdf", "Lognormal PDF", dsMass, mean3, sigma3)


  n_bkg_evt  = ROOT.RooRealVar(    "n_bkg_event",    "n_bkg_event", n_bkg)

  fracs = ROOT.RooRealVar(    "fracs",    "fracs",  0.5,     0.0,  1.0)
  fracs2 = ROOT.RooRealVar(    "fracs2",    "fracs2",  0.5,     0.0,  1.0)

  dg = ROOT.RooAddPdf( "dg","dg PDF",ROOT.RooArgList(gauss1,gauss2),fracs)
  gaussnormal = ROOT.RooAddPdf( "gaussnormal","gaussnormalPDF",ROOT.RooArgList(dg,lognormal),fracs2)
  #doublegauss = ROOT.RooExtendPdf("doublegauss","doublegauss PDF", dg,n_bkg_evt, "complete")

  #cheby = ROOT.RooChebychev ("cheby", "cheby", dsMass, ROOT.RooArgList(a1, a2, a3));
  #result2 = cheby.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #result2 = doublegauss.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #result2 = gauss1.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  result2 = gaussnormal.fitTo(hBkgRoo, ROOT.RooFit.Range("complete"), ROOT.RooFit.SumW2Error(True))
  #cheby.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("cheby"),             ROOT.RooFit.Components("cheby"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #doublegauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("doublegauss"),             ROOT.RooFit.Components("doublegauss"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )

  dsMass.setRange("left"     ,mlow3  ,mlow2)   
  dsMass.setRange("right"    ,mhigh2 ,mhigh3)        
  dsMass.setRange("signal"   ,mlow  ,mhigh)

  A = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"left"  ).getVal()
  C = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"right" ).getVal()
  B = n_bkg*gaussnormal.createIntegral(  dsMass,dsMass,"signal").getVal()
  import pdb
  print(f"integral values: A = {A}, B = {B}, C = {C}")

  #gauss1.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("gauss1"),             ROOT.RooFit.Components("gauss1"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #gaussnormal.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kGreen + 2),   ROOT.RooFit.Name("gaussnormal"),             ROOT.RooFit.Components("gaussnormal"),  ROOT.RooFit.FillColor(ROOT.kGreen + 2), )
  #create extended bkg hist
  #hist_bkg_ext = hist_bkg.Clone()
  #hist_bkg_ext.Scale(n_bkg)
  #dg_events = doublegauss.getVal()
  #print("dg events", dg_events)

  #create frame3
  frame3 = dsMass.frame(ROOT.RooFit.Title(var))
  frame3.GetXaxis().SetTitle("D_{s} mass")
  frame3.GetXaxis().SetTitleOffset(1.5)
  #frame3.GetYaxis().SetRangeUser(0, 200000)
  #max_y = hBkgR.GetMaximum() * 1.2  # Scale the max value slightly for better visualization
  frame3.GetYaxis().SetRangeUser(0, max_y)


  #first plot data and model, they define the normalization range 
  hBkgRoo.plotOn(frame3,ROOT.RooFit.Name("data2"),  ROOT.RooFit.SumW2Error(True)  )
  gaussnormal.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed -7 ),       ROOT.RooFit.Name("gaussnormal"), ROOT.RooFit.FillColor(ROOT.kRed  - 7), ROOT.RooFit.DrawOption("F") )
  hBkgRoo.plotOn(frame3,ROOT.RooFit.Name("data3"))
  ##get reduce chi2, 8 is the number of parameters
  #chi2  = round(frame3.chiSquare("pdf_total","hData_Roo", 2), 5)

  ##plot only RooFit method

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
  #
  ##plot subcomponents
  ##pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray),      ROOT.RooFit.FillColor(ROOT.kGray),      ROOT.RooFit.Name("pdf_bkg"),          ROOT.RooFit.Components("pdf_bkg"),         ROOT.RooFit.DrawOption("F"),  ROOT.RooFit.FillStyle(1001))
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.DrawOption("F") , ROOT.RooFit.FillStyle(3444)        )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kRed - 7 ),   ROOT.RooFit.Name("pdf_sig"),             ROOT.RooFit.Components("pdf_sig"), ROOT.RooFit.FillColor(ROOT.kRed - 7)        )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.FillStyle(3344)   )
  #pdf_total.plotOn(frame3,ROOT.RooFit.LineColor(ROOT.kGray + 2),   ROOT.RooFit.Name("pdf_bkg"),             ROOT.RooFit.Components("pdf_bkg"),  ROOT.RooFit.FillColor(ROOT.kGray + 2), )


  c2.cd()
  pad1.cd()
  frame3.Draw("SAME")

  c2.Modified()
  c2.Update()
  c2.SaveAs(outdir + "signflip_mass_bebi.pdf")




  print("==============> Saved!")
  return n_bkg, n_sig, [A,B,C]

