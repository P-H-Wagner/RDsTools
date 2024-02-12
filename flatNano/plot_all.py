
###############################
## Plotting script           ##
###############################

import ROOT

#input files
files = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/*"

#chain them
chain = ROOT.TChain("tree")
chain.Add(files)

#get branche names
names = [branch.GetName() for branch in chain.GetListOfBranches()]

#get number of entries
nEntries = chain.GetEntries();

#colors to multiVar
colors = [
    ROOT.kBlack,
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kOrange,
    ROOT.kMagenta,
    ROOT.kCyan,
    ROOT.kYellow,
    ROOT.kGray,
    ROOT.kViolet,
    ROOT.kTeal,
    ROOT.kSpring,
    ROOT.kAzure,
    ROOT.kPink,
    ROOT.kSpring + 6,
    ROOT.kTeal + 6,
    ROOT.kAzure + 6,
    ROOT.kPink + 6
]

#signal indices
sigs = [0,1,2,3]

#signal names
sigNames = {"0": r"D_{s} #mu","1": r"D_{s}* #mu", "2":r"D_{s} #tau", "3":r"D_{s}* #tau"}

#disable title and stats box
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

################################################################################################
#function which plots the desired variables into one canvas 

def multiVar(variables, bins, begin, end, name, norm = True, xLabel = "[GeV]", yLabel = "a.u."):

  """
  variables        = list of strings, the variables to plot
  begin, en d      = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.png)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  #define histograms
  for i,var in enumerate(variables):
    histos[var] = ROOT.TH1F(var,var,bins,begin,end)
    histos[var].SetLineColor(colors[i])
    histos[var].SetLineWidth(2)

  #fill all histograms

  print "Start filling " + name + " histogram; Events to fill: " + str(nEntries)

  for i in range(nEntries):
    if (i%20000 == 0):
      print " ==> Processing " + str(i) + "th Event"
    chain.GetEntry(i)
    for var in variables:
      histos[var].Fill(getattr(chain,var))

  #scale histograms
  if norm:
    for var in histos.keys():
      histos[var].Scale(1/histos[var].Integral())

  #get maximumm of y axis 
  yMax = max([histos[var].GetMaximum() for var in histos.keys()]) 

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  canvas.cd()
  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)

  for i,var in enumerate(variables):

    #multiVar
    if i == 0:
      histos[var].SetMaximum(yMax*1.2)
      histos[var].GetYaxis().SetTitle(yLabel)
      histos[var].GetXaxis().SetTitle(xLabel)
      histos[var].Draw("HIST")
    else:
      histos[var].Draw("HIST SAME")

    #legend
    legend.AddEntry(histos[var], var, "l")

  legend.Draw("SAME")

  #saving
  canvas.SaveAs("./plots/" + name + ".png")



################################################################################################
#function which plots the desired signals into one canvas 
def multiSig(var, bins, begin, end, name, signals = sigs, norm = True, xLabel = "[GeV]", yLabel = "a.u."):

  """
  variables        = list of strings, the variables to plot
  begin, en d      = floats, begin and end of range
  bins             = int, nr of bins
  name             = str, the name under which the plot should be saved (name.png)
  norm             = bool, tells if we normalize the histograms
  xLabel, yLabel   = strings, axis labels
  """

  #holds histograms
  histos = {}
  addresses = {}

  
  #define histograms
  for i,sig in enumerate(signals):
    histos[str(sig)] = ROOT.TH1F(str(sig),str(sig),bins,begin,end)
    histos[str(sig)].SetLineColor(colors[i])
    histos[str(sig)].SetLineWidth(2)

  #fill all histograms

  print "Start filling " + name + " histogram; Events to fill: " + str(nEntries)

  for i in range(nEntries):
    if (i%20000 == 0):
      print " ==> Processing " + str(i) + "th Event"
    chain.GetEntry(i)
    for sig in signals:
      if (sig == getattr(chain,"sig")):
        histos[str(sig)].Fill(getattr(chain,var))
  
  #scale histograms
  if norm:
    for sig in signals:
      histos[str(sig)].Scale(1/histos[str(sig)].Integral())

  #get maximumm of y axis 
  yMax = max([histos[str(sig)].GetMaximum() for sig in signals]) 

  #Draw histos
  canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
  canvas.cd()
  legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9); # Specify legend coordinates (x1, y1, x2, y2)
  legend.SetTextSize(0.03)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)

  for i,sig in enumerate(signals):

    #multiVar
    if i == 0:
      histos[str(sig)].SetMaximum(yMax*1.2)
      histos[str(sig)].GetYaxis().SetTitle(yLabel)
      histos[str(sig)].GetXaxis().SetTitle(xLabel)
      histos[str(sig)].Draw("HIST")
    else:
      histos[str(sig)].Draw("HIST SAME")

    #legend
    legend.AddEntry(histos[str(sig)], sigNames[str(sig)], "l")

  legend.Draw("SAME")

  #saving
  canvas.SaveAs("./plots/" + name + ".png")




###################################
## Call Plotting function        ##
###################################

multiSig("lxy_bs" ,50,0,3   ,"lxy_bs" ,xLabel = "[cm]")
multiSig("lxy_ds" ,50,0,1   ,"lxy_ds" ,xLabel = "[cm]")
multiSig("lxy_phi",50,0,1.8 ,"lxy_phi",xLabel = "[cm]")

multiSig("dxy_mu" ,50,0,0.3 ,"dxy_mu" ,xLabel = "[cm]")
multiSig("dxy_pi" ,50,0,0.6 ,"dxy_pi" ,xLabel = "[cm]")
multiSig("dxy_k1" ,50,0,0.4 ,"dxy_k1" ,xLabel = "[cm]")

multiSig("dz_mu"  ,50,0,0.3 ,"dz_mu" ,xLabel = "[cm]")
multiSig("dz_pi"  ,50,0,0.7 ,"dz_pi" ,xLabel = "[cm]")
multiSig("dz_k1"  ,50,0,0.5 ,"dz_k1" ,xLabel = "[cm]")

#bs momentum comparison
multiVar(["bs_pt_coll", "bs_pt_lhcb", "bs_pt_lhcb_alt", "bs_gen_pt", "bs_pt_reco_1", "bs_pt_reco_2"], 30,0.0,150.0, "bsComparison")

#hel angles

#collinear approx spoils this distribution as it ds mu system is back to back in bs rest frame --> dont plot it
multiVar(["cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen"], 30,-1,1, "cosMuWComparison",xLabel = "") 
multiVar(["cosPhiDsColl","cosPhiDsLhcb", "cosPhiDsLhcbAlt", "cosPhiDsReco1", "cosPhiDsReco2", "cosPhiDsGen"], 30,-1,1, "cosPhiDsComparison",xLabel = "") 
multiVar(["cosPiDsColl","cosPiDsLhcb", "cosPiDsLhcbAlt", "cosPiDsReco1", "cosPiDsReco2", "cosPiDsGen"], 30,-1,1, "cosPiDsComparison",xLabel = "") 
multiVar(["cosPlaneBsColl","cosPlaneBsLhcb", "cosPlaneBsLhcbAlt", "cosPlaneBsReco1", "cosPlaneBsReco2", "cosPlaneBsGen"], 30,-1,1, "cosPlaneBs",xLabel = "")
multiVar(["cosPlaneDsColl","cosPlaneDsLhcb", "cosPlaneDsLhcbAlt", "cosPlaneDsReco1", "cosPlaneDsReco2", "cosPlaneDsGen"], 30,-1,1, "cosPlaneDs",xLabel = "")

multiVar(["cosPiK1","cosPiK1Gen"],30,-1,1,"cosPiK1", xLabel = "")

#compare prefit, postfit and gen momenta
multiVar(["k1_pt",    "k1_refitted_pt",    "k1_gen_pt"],                     30,0,30, "pre_vs_postfit_k1")
multiVar(["k2_pt",    "k2_refitted_pt",    "k2_gen_pt"],                     30,0,30, "pre_vs_postfit_k2")
multiVar(["pi_pt",    "pi_refitted_pt",    "pi_gen_pt"],                     30,0,30, "pre_vs_postfit_pi")
multiVar(["mu_pt",    "mu_refitted_pt",    "mu_gen_pt"],                     30,0,30, "pre_vs_postfit_mu")
multiVar(["kk_pt",    "phi_fitted_pt",     "phi_refitted_pt", "phi_gen_pt"], 30,0,30, "pre_vs_postfit_phi")
multiVar(["phiPi_pt", "ds_fitted_pt",      "ds_refitted_pt", "ds_gen_pt"],   30,0,30, "pre_vs_postfit_ds")

