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

# PERFROM SIDEBAND RECONSTRUCTRION
# FIND: A,B,C,S


#############################################
#                                           #
#                                           #
#  *                                        #
#   *            *                          #
#     *        *   *                        #
#        *    *     *                       #
#           *    S   *                      #
#           |  *      *                     #
#           |      *   *                    #
#           |          | *                  #
#           |          |     *              #
#           |          |         *          #
#    A      |     B    |     C       *   *  #
#           |          |                    #
#############################################  Ds mass



############## FIT TO MC ################
 
# Binning:
bins  = 21
start = 1.91#1.91
stop  = 2.028#2.028

# Range for mean scale
meanScaleStart = 1.001 # dont start at 1
meanScaleMin   = 0.9
meanScaleMax   = 1.1

# Range for CB sigma
cbStart = 0.008
cbMin   = 0.005
cbMax   = 0.012

# Range for 1st gauss sigma scale
sigma1Start = 1.1
sigma1Min   = 0.05
sigma1Max   = 3

# Range for 2nd gauss sigma scale
sigma2Start = 0.99
sigma2Min   = 0.5
sigma2Max   = 100

# Range for alpha of CB
alphaStart  = 1.1
alphaMin    = 0.01
alphaMax    = 5

# Range for n of CB
nStart      = 2.1
nMin        = 0.5
nMax        = 7

# Range for 1st gauss weight
fracs1Start = 0.2
fracs1Min   = 0.1
fracs1Max   = 1

# Range for CB weight
fracs2Start = 0.9
fracs2Min   = 0.1
fracs2Max   = 1

# Nr of parameters in double gauss + cb (for red chi2)
nParams     = 8

########### FIT TO DATA ############

# Range of exp contant
aStart      = -2
aMin        = -4
aMax        = 0

# Range of sigma for data gauss
sigma3Start = 0.01
sigma3Min   = 0.005
sigma3Max   = 0.20

# Range of scale for data gauss sigma
scaleStart  = 2
scaleMin    = 1.1
scaleMax    = 4

# Range for data gauss weight

fracs3Start = 0.5 
fracs3Min   = 0.0
fracs3Max   = 1.0

# Range for data double gauss weight

fracs4Start = 0.8
fracs4Min   = 0.0
fracs4Max   = 1.0


# Selection
selec = ' & '.join([
f'dsMu_m < {bsMass_}',
'k1_charge*k2_charge <0',
'mu_charge*pi_charge < 0',
'mu_pt > 8',
'k1_pt > 1',
'k2_pt > 1',
'pi_pt > 1',
'lxy_ds < 1',
'mu_id_medium == 1',
'mu_rel_iso_03 < 0.3'
#'tv_prob > 0.1',
#'fv_prob > 0.1'
])


## Cuts from MA:
# 'mu_pt>8.0',
# 'k1_pt>1.0',
# 'k2_pt>1.0',
# 'pi_pt>1.0',
# 'HLT_Mu7_IP4',
# 'phi_vtx_prob>0.1',
# 'mu_id_medium',
# 'mu_rel_iso<0.3',
# 'lxy_ds < 1.0',
# 'k1_charge*k2_charge <0'])


selection     = selec + "& gen_sig == 0"
selectionData = selec 

#plotting settings
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True) #important for error calc.!
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadLeftMargin(0.15)


outdir = "./sb_fits/"


#importing trees

def getRdf(dateTime):

  files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{dateTime}/*" #test
  print(files)

  chain = ROOT.TChain("tree")
  chain.Add(files)
  rdf = ROOT.RDataFrame(chain)

  return (chain,rdf)

def defineRanges(var,nBins, start, stop):

  binEdges = np.linspace(start, stop, nBins+1)
  print(binEdges) 
  for i in range(nBins):
    #i starts at 0
    var.setRange("i+1", binEdges[i], binEdges[i+1])
  return var, binEdges

def getSigma(rdf, var, sel, bins = bins, start = start, stop = stop):

  print(" ======> Prepare variables and fit")
  #path to save
  os.system(f"mkdir -p ./plots/")

  model = ROOT.RDF.TH1DModel( var, "", bins, start, stop)

  #create histo
  h = rdf.Filter(sel).Histo1D(model,var)

  # adjust histo style
  h.GetXaxis().SetTitle(r"D_{s} mass [GeV]")
  h.GetXaxis().SetTitleOffset(1.3)
  h.GetYaxis().SetTitleSize(0.8*h.GetYaxis().GetTitleSize())
  h.GetXaxis().SetTitleSize(0.8*h.GetXaxis().GetTitleSize())
  h.GetYaxis().SetTitle('Events')
  h.SetLineColor(1)
  h.SetLineWidth(2)

  # and create ptr to this histo ( for the cosntructor of hRoo)
  hPtr = h.GetPtr()

  h.Scale(0.999) #apply?

  # Fitting variable
  dsMass = ROOT.RooRealVar(var,r"D_{s} mass", start, stop)

  # Define a RooDataHist from TH1D
  hRoo   = ROOT.RooDataHist("hRoo", "hRoo", ROOT.RooArgList(dsMass), hPtr)


  # Mean and scaled mean
  mean       = ROOT.RooRealVar(    "mean",       "mean",              dsMass_                                    ) 
  meanScale  = ROOT.RooRealVar(    "mean_scale", "mean_scale",        meanScaleStart, meanScaleMin, meanScaleMax ) 
  scaledMean = ROOT.RooFormulaVar( "scaledMean", "mean * mean_scale", ROOT.RooArgList(mean, meanScale)           ) 

  # Sigma for cb
  sigmaCb    = ROOT.RooRealVar(    "sigmaCb",    "sigmaCb",           cbStart,        cbMin,        cbMax        ) 

  # Sigma for first gaussian
  scale1     = ROOT.RooRealVar(    "scale1",     "scale1",            sigma1Start,    sigma1Min,    sigma1Max    )
  sigma1     = ROOT.RooFormulaVar( "sigma1",     "sigmaCb * scale1",  ROOT.RooArgList(sigmaCb,scale1)            )

  # Sigma for second gaussian
  scale2     = ROOT.RooRealVar(    "scale2",     "scale2",            sigma2Start,    sigma2Min,    sigma2Max    )
  sigma2     = ROOT.RooFormulaVar( "sigma2",     "sigmaCb * scale2",  ROOT.RooArgList(sigmaCb,scale2)            )


  #alpha and n for cb
  alpha      = ROOT.RooRealVar(    "alpha",      "alpha",             alphaStart,     alphaMin,     alphaMax     )
  n          = ROOT.RooRealVar(    "n",          "n",                 nStart,         nMin,         nMax         )

  # Define functions

  # define CB
  cb          = ROOT.RooCBShape("cb",     "cb",     dsMass,      scaledMean, sigmaCb, alpha, n)

  # define gaussians for double gauss
  gauss1      = ROOT.RooGaussian(   "gauss1", "gauss1", dsMass,      scaledMean, sigma1           )
  gauss2      = ROOT.RooGaussian(   "gauss2", "gauss2", dsMass,      scaledMean, sigma2           )

  # relative fractions between CB and double gauss
  fracs       = ROOT.RooRealVar(    "fracs",  "fracs",  fracs1Start, fracs1Min,  fracs2Max        ) #weight of gauss 1
  fracs2      = ROOT.RooRealVar(    "fracs2", "fracs2", fracs2Start, fracs2Min,  fracs2Max        ) #weight of CB
  
  # Define double gauss and CB+double gauss
  #doubleGauss = ROOT.RooAddPdf(     "doubleGauss", "doubleGauss", ROOT.RooArgList(gauss1,gauss2),  ROOT.RooArgList(fracs) )
  cbGauss    = ROOT.RooAddPdf(      "cbGauss",     "cbGauss",     ROOT.RooArgList(cb,gauss1), ROOT.RooArgList(fracs2))
  
  #perform the fit
  dsMass.setRange("dsMassRange", start, stop)
  result = cbGauss.fitTo(hRoo, ROOT.RooFit.Range("dsMassRange"), ROOT.RooFit.SumW2Error(True))
  
  #get results of parameters
  alphaCB = round(alpha.getVal()      ,5)
  nCb     = round(n.getVal()          ,5)
  s1      = round(sigmaCb.getVal()    ,5) #-->the sigma for the SB
  f1      = round(meanScale.getVal() ,5)
  rel2    = round(fracs2.getVal()     ,2)
  
  print(f" =======> simga from cb = {s1}, weight of cb = {rel2}")

  #plotting
  print(f" =======> Start plotting")

  #create frame
  frame = dsMass.frame(ROOT.RooFit.Title(var))
  frame.GetXaxis().SetTitle("D_{s} mass")
  frame.GetXaxis().SetTitleOffset(1.5)

  #first plot data and model, they define the normalization range 
  hRoo.plotOn(frame,ROOT.RooFit.Name("data"))
  cbGauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("cbGauss")) #,ROOT.RooFit.Components("cbGauss"))

  #get reduce chi2, 8 is the number of parameters
  chi2  = round(frame.chiSquare("cbGauss","data", nParams), 5)

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
  #cbGauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kCyan +1), ROOT.RooFit.Name("doubleGauss"), ROOT.RooFit.Components("doubleGauss"))
  cbGauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kCyan +1), ROOT.RooFit.Name("gauss1"), ROOT.RooFit.Components("gauss1"))
  cbGauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue ),   ROOT.RooFit.Name("cb"),          ROOT.RooFit.Components("cb")         )
  cbGauss.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("cbGauss")) #,ROOT.RooFit.Components("cbGauss"))

  #plot data again, s.t. it is visible on top
  hRoo.plotOn(frame,ROOT.RooFit.Name("data"))
  
  #get pulls and add them to frame2
  hpull = frame.pullHist("data","cbGauss")
  frame2 = dsMass.frame(ROOT.RooFit.Title(""))
  frame2.addPlotable(hpull,"P")
  
  c2.cd()
  pad1.cd()
  frame.Draw()
  
  leg2 = ROOT.TLegend(0.6,0.7,0.8,0.88);
  leg2.SetFillColor(ROOT.kWhite);
  leg2.SetLineColor(ROOT.kWhite);
  leg2.SetTextSize(0.030)
  leg2.AddEntry("data","MC Data","EP")
  leg2.AddEntry("cb","Crystal Ball (CB)","L")
  leg2.AddEntry("gauss1","Gaussian (G)","L")
  leg2.AddEntry("cbgauss","CB + G","L")
  leg2.Draw("SAME")
  
  text = ROOT.TPaveText(.18,.7,.4,.88,"brNDC");
  text.SetFillColor(0)
  text.SetFillStyle(0)
  text.SetTextAlign(11)
  leg2.SetTextSize(0.028)
  text.AddText("#sigma_{CB} = " + f"{s1} GeV")
  text.AddText("#mu_{m} = " + f"{f1} ")
  text.AddText("#chi^{2} = " +f"{round(chi2,2)}")
  text.SetBorderSize(0)
  text.Draw("SAME")
  
  CMS_lumi(pad1, 4, 0, cmsText = 'CMS', extraText = '  Private work', lumi_13TeV = '')
  
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
  	event = hpull.Eval(h.GetBinCenter(i+1))
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
  c2.SaveAs(outdir + "RooFit_only.pdf")


  return s1, h


######################################################################################

def getABCS(filename, rdf, sel, var, sigma, hMc, bins = bins, start = start, stop = stop, binsFake = 21, nSig = nSignalRegion, nSb = nSidebands, width = sbWidth):

  """
  - sigma obtained from fit to MC"  
  - defines width of signal region (in terms of sigmas)
  - defines width of sb region (in terms of sigmas)
  """

  print("========> Get ABCS")


  #signal region
  mlow   = dsMass_ - nSig*sigma
  mhigh  = dsMass_ + nSig*sigma

  #sideband start
  mlow2  = dsMass_ - nSb*sigma
  mhigh2 = dsMass_ + nSb*sigma

  #sideband stops
  mlow3  = mlow2  - width*sigma
  mhigh3 = mhigh2 + width*sigma 

  #create histo

  model = ROOT.RDF.TH1DModel( var, "", bins, start, stop)
  h = rdf.Filter(sel).Histo1D(model,var)

  # adjust histo style
  h.GetXaxis().SetTitle(r"D_{s} mass [GeV]")
  h.GetXaxis().SetTitleOffset(1.3)
  h.GetYaxis().SetTitleSize(0.8*h.GetYaxis().GetTitleSize())
  h.GetXaxis().SetTitleSize(0.8*h.GetXaxis().GetTitleSize())
  h.GetYaxis().SetTitle('Events')
  h.SetMarkerColor(ROOT.kBlack)
  h.SetFillColor(ROOT.kBlack)
  h.SetLineColor(ROOT.kBlack)
  h.SetMarkerStyle(8)	

  # and a pointer to the histogram for the RooDataHist  
  hPtr = h.GetPtr()

  #total number of data events
  #nData = h.GetEntries()
  nData = h.Integral()
  print("calculating fake mass with nData:", nData) 

 
  #observable is again the Ds mass
  combMass = ROOT.RooRealVar("combMass",r"D_{s} mass",start,stop)
  
  #observable ranges
  combMass.setRange("left"     ,mlow3  ,mlow2)   
  combMass.setRange("right"    ,mhigh2 ,mhigh3)        
  combMass.setRange("signal"   ,mlow  ,mhigh)
  combMass.setRange("complete" ,start ,stop)   
  combMass.setRange("full"     ,mlow3 ,mhigh3)
  
  #Create edges for the fake mass histogram, s.t. we can plot a MC Ds mass distribution
  binEdges = np.linspace(start, stop, binsFake+1)
  #print(binEdges) 
  # Defien ranges for every bin, this does not work within a for loop -.-
  #for i in range(binsFake):
  #  #i starts at 0
  #  combMass.setRange("i+1", binEdges[i], binEdges[i+1])


  combMass.setRange("1",  binEdges[0],  binEdges[1])
  combMass.setRange("2",  binEdges[1],  binEdges[2])
  combMass.setRange("3",  binEdges[2],  binEdges[3])
  combMass.setRange("4",  binEdges[3],  binEdges[4])
  combMass.setRange("5",  binEdges[4],  binEdges[5])
  combMass.setRange("6",  binEdges[5],  binEdges[6])
  combMass.setRange("7",  binEdges[6],  binEdges[7])
  combMass.setRange("8",  binEdges[7],  binEdges[8])
  combMass.setRange("9",  binEdges[8],  binEdges[9])
  combMass.setRange("10", binEdges[9],  binEdges[10])
  combMass.setRange("11", binEdges[10], binEdges[11])
  combMass.setRange("12", binEdges[11], binEdges[12])
  combMass.setRange("13", binEdges[12], binEdges[13])
  combMass.setRange("14", binEdges[13], binEdges[14])
  combMass.setRange("15", binEdges[14], binEdges[15])
  combMass.setRange("16", binEdges[15], binEdges[16])
  combMass.setRange("17", binEdges[16], binEdges[17])
  combMass.setRange("18", binEdges[17], binEdges[18])
  combMass.setRange("19", binEdges[18], binEdges[19])
  combMass.setRange("20", binEdges[19], binEdges[20])
  combMass.setRange("21", binEdges[20], binEdges[21])

  #create a Roo histo from our TH1D histo
  hRoo   =  ROOT.RooDataHist("hRoo","hRoo",ROOT.RooArgList(combMass),hPtr)
  
  #make it a RooRealVar but it is fixed
  nDataVar       = ROOT.RooRealVar(    "nDataVar",       "nDataVar",             nData) 
  #parameters for exp
  a      = ROOT.RooRealVar(      "a",     "a",     aStart, aMin, aMax)

  #parameters for double gauss
  sigma3 = ROOT.RooRealVar(      "sigma3","sigma3", sigma3Start, sigma3Min, sigma3Min)
  scale  = ROOT.RooRealVar(      "scale","scale", scaleStart, scaleMin, scaleMax)
  sigma4 = ROOT.RooFormulaVar(   "sigma4","sigma3 * scale",ROOT.RooArgList(sigma3,scale))
  mean   = ROOT.RooRealVar(      "mean","mean",dsMass_)

  #define gaussians fro double gauss
  gauss3 = ROOT.RooGaussian(     "gauss3", "gaussian PDF1", combMass, mean, sigma3)
  gauss4 = ROOT.RooGaussian(     "gauss4", "gaussian PDF2", combMass, mean, sigma4)

  #define relative weights 
  fracs3 = ROOT.RooRealVar(      "fracs3",  "fracs3",         fracs3Start, fracs3Min, fracs3Max)   #weight for gauss 3
  fracs4 = ROOT.RooRealVar(      "fracs4", "fracs4",        fracs4Start, fracs4Min, fracs4Max) #weight of double gauss
  
  #define double gauss, exp, and double gauss + exp
  doubleGauss = ROOT.RooAddPdf( "doublegauss","doublegauss PDF",ROOT.RooArgList(gauss3,gauss4),fracs3)

  #extend the gauss!
  expSig        = ROOT.RooRealVar(      "expSig","expSig", nData * 0.07, 0, nData) # nr of signals
  doubleGaussEx = ROOT.RooExtendPdf("ext", "ext", doubleGauss,expSig , "complete") # extended sig pdf

  exp = ROOT.RooExponential(     "exp","exp PDF",combMass,a)
  #extend the exp!
  expBkg        = ROOT.RooFormulaVar(      "expBkg","nDataVar - expSig",ROOT.RooArgList(nDataVar,expSig) ) # nr of bkg
  expEx = ROOT.RooExtendPdf("ext2", "ext2",exp,expBkg , "complete") # extended bkg pdf

  model = ROOT.RooAddPdf(        "model","model PDF",ROOT.RooArgList(doubleGauss,exp),fracs4)
  #EXtend the model!
  modelEx = ROOT.RooAddPdf(        "modelEx","model PDFEx",ROOT.RooArgList(doubleGaussEx,expEx),ROOT.RooArgList(expSig, expBkg))
  
 
  #perform fit
  #result = model.fitTo(hRoo)     
  resultEx = modelEx.fitTo(hRoo)     
  print("extended fitted nsig = ", expSig.getVal())

  #extract number of signals and bkg
  #print("fracs4: ",fracs4.getVal())
  #nSig = nData * fracs4.getVal()
  #nBkg = nData - nSig
  nSig = expSig.getVal()
  nBkg = expBkg.getVal()  

  print("test: nSig + nBkg = ", nSig + nBkg, "= nData = ", nData)

  #define frame to plot on
  frame_comb = combMass.frame(ROOT.RooFit.Title(""))
  frame_comb.GetXaxis().SetTitle(r"D_{s} mass")
  frame_comb.GetXaxis().SetTitleOffset(1.5)
  frame_comb.GetYaxis().SetNdivisions(10)
  frame_comb.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.4)
  
  #plotting
  c_comb = ROOT.TCanvas("RooComb","RooComb",700,700)
  pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.2, 0.89, 0.99)
  pad1.SetLeftMargin(0.15)
  pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.89, 0.2)
  pad2.SetLeftMargin(0.15)
  pad3 = ROOT.TPad("pad3", "pad3", 0.82, 0.01, 0.99, 0.2)
  pad3.SetLeftMargin(0.15)
  
  pad1.Draw()
  pad2.Draw()
  pad3.Draw()
  
  
  #print("first fake before drawing:", nBkg  * exp.createIntegral(  combMass, combMass, "1"  ).getVal())
  

  hRoo.plotOn(frame_comb,ROOT.RooFit.Name("data"))#,ROOT.RooFit.Range("complete")
  modelEx.plotOn(frame_comb,ROOT.RooFit.LineColor(ROOT.kRed -7),ROOT.RooFit.FillColor(ROOT.kRed -7),ROOT.RooFit.DrawOption('F'),ROOT.RooFit.LineStyle(10),ROOT.RooFit.Components("modelEx"),ROOT.RooFit.Name("modelEx"))#,ROOT.RooFit.Range("complete"))
  modelEx.plotOn(frame_comb,ROOT.RooFit.LineColor(ROOT.kGray ),ROOT.RooFit.FillColor(ROOT.kGray ),ROOT.RooFit.DrawOption('F'),ROOT.RooFit.LineStyle(0),ROOT.RooFit.Components("ext2"),ROOT.RooFit.Name("ext2")) #,ROOT.RooFit.Range("complete"))
  #model.plotOn(frame_comb,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.LineStyle(10),ROOT.RooFit.Components("doublegauss"),ROOT.RooFit.Name("double gauss"))#,ROOT.RooFit.Range("complete"))
  hRoo.plotOn(frame_comb,ROOT.RooFit.Name("data"))#,ROOT.RooFit.Range("complete")

  # we have to call createIntegral at this part of the code, otherwise it will crash!!
  
  #get normalization constants A,B,C,S

  #sidebands
  #A = nBkg  * exp.createIntegral(  combMass,combMass,"left"  ).getVal()
  #C = nBkg  * exp.createIntegral(  combMass,combMass,"right" ).getVal()
  A = nBkg*expEx.createIntegral(  combMass,combMass,"left"  ).getVal()
  C = nBkg*expEx.createIntegral(  combMass,combMass,"right" ).getVal()

  #signal and backgroudn entries
  #B = nBkg  * exp.createIntegral(  combMass,combMass,"signal").getVal()
  #S = nData * model.createIntegral(combMass,combMass,"signal").getVal() - nBkg*exp.createIntegral(combMass,combMass,"signal").getVal()
  B = nBkg*expEx.createIntegral(  combMass,combMass,"signal").getVal()
  S = nData*modelEx.createIntegral(combMass,combMass,"signal").getVal() - nBkg*expEx.createIntegral(combMass,combMass,"signal").getVal()

  print(f"integral values: A = {A}, B = {B}, C = {C}, S = {S}")
  print(f"From the fit we expect {A+B+C} comb. bkg. events int the complete region and {B} events in the signal region")
  # get bin content for the fake mass
  bin_content = []
  
  for i in range(binsFake):
    nFakes = int(nBkg* expEx.createIntegral(combMass,combMass,f"{i+1}").getVal())
    nSigFakes = int(nSig* doubleGaussEx.createIntegral(combMass,combMass,f"{i+1}").getVal())
    print("data is:", h.GetBinContent(i+1))
    print("data edge is:", h.GetXaxis().GetBinLowEdge(i+1))
    print("sig is:", nSigFakes)
    print("fakes are:", nFakes)
    bin_content.append(nFakes)
  
  #continue plotting....
  
  #get pulls
  hpull = frame_comb.pullHist("data","modelEx")
  
  frame_comb2 = combMass.frame(ROOT.RooFit.Title(""))
  frame_comb2.addPlotable(hpull,"P")
  frame_comb2.GetYaxis().SetTitle(r"#Delta / #sigma")
  frame_comb2.GetXaxis().SetTickLength(0.1)
  frame_comb2.GetXaxis().SetTitleSize(3.* frame_comb2.GetXaxis().GetTitleSize())
  frame_comb2.GetYaxis().SetTitleSize(4.* frame_comb2.GetYaxis().GetTitleSize())
  frame_comb2.GetXaxis().SetLabelSize(3.* frame_comb2.GetXaxis().GetLabelSize())
  frame_comb2.GetYaxis().SetLabelSize(3.* frame_comb2.GetYaxis().GetLabelSize())
  frame_comb2.GetYaxis().SetTitleOffset( 0.34 )
  frame_comb2.GetYaxis().CenterTitle()
  ratio_bins = 11
  frame_comb2.GetYaxis().SetNdivisions(ratio_bins)
  ratio_start = -4
  ratio_stop = 4
  frame_comb2.GetYaxis().SetRangeUser(ratio_start,ratio_stop)
  
  #get chi2, we have 5 parameters
  chi2 =round(frame_comb.chiSquare("modelEx","data",5),5)
  
  c_comb.cd()
  pad1.cd()
  
  #extract exp as TF1 
  #expTF1 = frame_comb.findObject("exp")
  expTF1 = frame_comb.findObject("ext2")

  #print(exp.createIntegral(  combMass,combMass).getVal())
  #width = binEdges[1]-binEdges[0]
  #exp2  = exp.asTF(ROOT.RooArgList(combMass), ROOT.RooArgList(a), ROOT.RooArgList(combMass))
  #mod2 = model.asTF(ROOT.RooArgList(combMass), ROOT.RooArgList(a), ROOT.RooArgList(combMass))

  #print("nBkg", nBkg)
  #print("bin1 with TF1", nBkg*exp2.Integral(binEdges[0],  binEdges[1]))
 
  #create simple lienar TF1 (approximation to exp)
  binwidth = h.GetBinWidth(1)
  
  #to visualize sidebands ;)
  x1 = mlow2
  x2 = mhigh2
  y1 = expTF1.Eval(mlow2)
  y2 = expTF1.Eval(mhigh2)
  slope = (y2-y1)/(x2-x1) 
  b = y1 - slope*x1
  dummy1 = ROOT.TF1("dummy1","[0] + x*[1] ",mlow2,mlow3)
  dummy2 = ROOT.TF1("dummy2","[0] + x*[1] ",mhigh2,mhigh3)
  dummy1.SetParameter(0,0.999*b)
  dummy1.SetParameter(1,slope)
  dummy2.SetParameter(0,0.999*b)
  dummy2.SetParameter(1,slope)
  dummy1.SetFillStyle(3244)
  dummy2.SetFillStyle(3244)
  dummy1.SetFillColor(ROOT.kBlue)
  dummy2.SetFillColor(ROOT.kBlue)
  dummy1.SetLineColor(ROOT.kBlue)
  dummy2.SetLineColor(ROOT.kBlue)
  frame_comb.Draw()
  dummy1.Draw("SAME F")
  dummy2.Draw("SAME F")
  h.Draw("SAME EP")
 
  #legend
  legcomb = ROOT.TLegend(.2,.8,.6,.88)
  legcomb.SetBorderSize(0)
  legcomb.SetFillColor(0)
  legcomb.SetFillStyle(0)
  legcomb.SetTextFont(42)
  legcomb.SetNColumns(2)
  legcomb.SetTextSize(0.035)
  legcomb.AddEntry("model" ,'Signal','F')
  legcomb.AddEntry("data" ,'Data' ,'LEP' )
  legcomb.AddEntry("exp" ,'Comb. Bkg.','F')
  legcomb.AddEntry("dummy1","Sidebands", "F")
  legcomb.Draw("SAME")
  
  #text
  textcomb = ROOT.TPaveText(.65,.82,.75,.88,"brNDC");
  textcomb.SetTextFont(43)
  textcomb.SetFillColor(0)
  textcomb.SetFillStyle(0)
  textcomb.SetBorderSize(0)
  textcomb.SetTextSize(16)
  textcomb.AddText("#chi^{2}_{red}" + f" = {round(chi2,2)}")
  textcomb.Draw("SAME")
  
  CMS_lumi(pad1, 4, 0, cmsText = '     CMS', extraText = '       Preliminary', lumi_13TeV = 'X fb^{-1}')
  
  pad2.cd()
  
  frame_comb2.Draw()
  
  #add lines at +-1 at resiudal plot
  l1 = ROOT.TLine(start,1,stop,1)
  l2 = ROOT.TLine(start,-1,stop,-1)
  l1.SetLineStyle(2)
  l2.SetLineStyle(2)
  l1.Draw()
  l2.Draw()
  
  pad3.cd()
  
  #prepare hsitogram for residual distribution
  hgauss = ROOT.TH1D("gauss","gauss",ratio_bins,ratio_start,ratio_stop)
  for i in range(bins):
  	#fill hgauss
  	event = hpull.Eval(hMc.GetBinCenter(i+1))
  	a = hgauss.Fill(event)
  hgauss.GetYaxis().SetTickSize(0)
  hgauss.GetXaxis().SetTickSize(0.1)
  hgauss.GetXaxis().SetNdivisions(ratio_bins)
  hgauss.GetYaxis().SetLabelSize(0)
  hgauss.GetXaxis().SetLabelSize(3.* hgauss.GetXaxis().GetLabelSize())
  hgauss.SetFillColor(ROOT.kBlack)
  hgauss.SetLineColor(ROOT.kBlack)
  hgauss.Draw("hbar")
  
  c_comb.Modified()
  c_comb.Update()
  c_comb.SaveAs(outdir + "expfit.pdf")
  ########################################################################################
  print("========> build fake mass histogram")
  
  #lets create a fake array
  fake_mass = []
  
  for i,events in enumerate(bin_content):
  
  	#this mass lies within the bin
  	mass_average = (binEdges[i] + binEdges[i+1])/2 
  	print(mass_average) 
  	for j in range(events):
  		#lets append this mass events-times
  		fake_mass.append(mass_average)
  
  #consistency check:
  check = len(fake_mass) - sum(np.array(bin_content))
  print(f"we are consisten if check = 0, and check is: {check}")
  
  #save bin content and fakemass histogram as csv
  np.savetxt("mass_bincontent.csv",bin_content,delimiter = ",")
  np.savetxt("fakemass.csv", fake_mass, delimiter=",")
  

  return A, B, C, S

"""

# this is for a test 
test = "26_04_2024_16_28_22"
ch, df = getRdf(test)
s, h  = getSigma(df, "phiPi_m", selection)
                                 #old uncosntr          #new constr
dataTest = "26_04_2024_18_08_15" #"14_04_2024_22_09_50" #"19_04_2024_17_22_22"#"10_04_2024_00_32_58" bad stats!
chData,dfData = getRdf(dataTest)
A, B, C, S    = getABCS( dataTest, dfData, selectionData, "phiPi_m", s, h, bins = 15) 

"""
