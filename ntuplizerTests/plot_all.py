
###############################
## Plotting script           ##
###############################

import ROOT

#input files
files = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/10_02_2024_19_14_59/*"

#chain them
chain = ROOT.TChain("Events")
chain.Add(files)

#get branche names
names = [branch.GetName() for branch in chain.GetListOfBranches()]

#get number of entries
nEntries = chain.GetEntries();

#colors to draw
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

#disable title and stats box
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

##############################################################
#function which plots the desired variables into one canvas ##
##############################################################


def draw(variables, bins, begin, end, name, norm = True, xLabel = "[GeV]", yLabel = "a.u."):

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
      histos[var].Fill(getattr(chain,var)[0])

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

    #draw
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
  canvas.SaveAs("./" + name + ".png")


###################################
## Call Plotting function        ##
###################################

#bs momentum comparison
draw(["bs_pt_coll", "bs_pt_lhcb", "bs_pt_lhcb_alt", "bs_gen_pt", "bs_pt_reco_1", "bs_pt_reco_2"], 30,0.0,150.0, "bsComparison")
#hel angles
draw(["cosMuWColl", "cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen"], 30,-1,1, "cosMuWComparison",xLabel = "")
draw(["cosPhiDs", "cosPhiDsGen"],30,-1,1,"cosPhiDs")

#compare prefit, postfit and gen momenta
draw(["k1_pt",    "k1_refitted_pt",    "k1_gen_pt"],                     30,0,30, "pre_vs_postfit_k1")
draw(["k2_pt",    "k2_refitted_pt",    "k2_gen_pt"],                     30,0,30, "pre_vs_postfit_k2")
draw(["pi_pt",    "pi_refitted_pt",    "pi_gen_pt"],                     30,0,30, "pre_vs_postfit_pi")
draw(["mu_pt",    "mu_refitted_pt",    "mu_gen_pt"],                     30,0,30, "pre_vs_postfit_mu")
draw(["kk_pt",    "phi_fitted_pt",     "phi_refitted_pt", "phi_gen_pt"], 30,0,30, "pre_vs_postfit_phi")
draw(["phiPi_pt", "ds_fitted_pt",      "ds_refitted_pt", "ds_gen_pt"],   30,0,30, "pre_vs_postfit_ds")


