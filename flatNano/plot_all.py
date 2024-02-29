
###############################
## Plotting script           ##
###############################

import ROOT
import math
import numpy as np

import argparse
parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help') 

parser.add_argument('filename')   
args = parser.parse_args()

#input files
files = f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/{args.filename}/*"

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

labelVar={"q2": "Q^{2} ", 
        "m2_miss": r"m^{2}_#mathrm{miss}",
        "cosMuW": r"cos(#mu W)",
        "cosPhiDs": r"cos(#phi D_{s})",
        "cosPiDs": r"cos(#pi D_{s})",
}

labelMethod={
        "coll": r"Collinear",
        #"lhcb": r"LHCb z",
        "lhcb_alt": r"LHCb xyz",
        "reco_1": r"Math Sol. 1", 
        "reco_2": r"Math Sol. 2",
        "gen"   : r"Gen",
        "Coll": r"Collinear",
        #"lhcb": r"LHCb z",
        "LhcbAlt": r"LHCb xyz",
        "Reco1": r"Math Sol. 1", 
        "Reco2": r"Math Sol. 2",
        "Gen"   : r"Gen",

}

#signal indices
sigs = [0,1,5,6]

#signal names
sigNames = {"0": r"D_{s} #mu","1": r"D_{s}* #mu", "2":r"D_{s} #tau", "3":r"D_{s}* #tau"}

#disable title and stats box
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

################################################################################################
#function which plots the desired variables into one canvas 

def multiVar(variables, bins, begin, end, name, norm = True, xLabel = "[GeV]", yLabel = "a.u.",color = colors):

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
    histos[var].SetLineColor(color[i])
    histos[var].SetLineWidth(2)

    histos["best"] = ROOT.TH1F("best","best",bins,begin,end)
    histos["best"].SetLineColor(ROOT.kViolet - 9)
    histos["best"].SetLineWidth(2)

  #fill all histograms

  print("Start filling " + name + " histogram; Events to fill: " + str(nEntries))

  for i in range(100):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)

    reco1Val = np.nan
    reco2Val = np.nan
    genVal   = np.nan
 
    for var in variables:
      if (getattr(chain,"sig") == 0):
        #if not (math.isnan(getattr(chain,"bs_pt_reco_1"))):
          histos[var].Fill(getattr(chain,var))

          #Get the recos and gen methods
          if ("reco_1" in var) or ("Reco1" in var): reco1Val = getattr(chain,var)
          if ("reco_2" in var) or ("Reco2" in var): reco2Val = getattr(chain,var)
          if ("gen" in var):                        genVal   = getattr(chain,var)

    

    # get the better solution (closer to gen)
    recos = [reco1Val, reco2Val]
    diffs = [abs(reco1Val - genVal), abs(reco2Val - genVal)]
   
    # Fillt histo with the better solution 
    histos["best"].Fill(recos[np.argmin(diffs)])


  #scale histograms
  if norm:
    for var in histos.keys():
      histos[var].Scale(1/histos[var].Integral())
    histos["best"].Scale(1/histos["best"].Integral()) ##CHANGE 

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


    # legend
    varLegend = r""

    #loop over methods
    for iName in labelMethod.keys():
      if iName in var: varLegend += labelMethod[iName]
    if ("lhcb" in var and "lhcb_alt" not in var) or ("Lhcb" in var and "LhcbAlt" not in var): varLegend += "LHCb z" #special bc lhcb also in lhcb_alt

    legend.AddEntry(histos[var], varLegend, "l")

   
  histos["best"].Draw("HIST SAME") ##CHANGE
  legend.AddEntry(histos["best"], "Math Sol. Best", "l")


  legend.Draw("SAME")

  #saving
  canvas.SaveAs("./plots/" + name + ".png")



################################################################################################
#function which plots the desired signals into one canvas 
def multiSig(var, bins, begin, end, name, signals = sigs, norm = True, xLabel = r" Q^{2} [GeV]", yLabel = "a.u."):

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

  print("Start filling " + name + " histogram; Events to fill: " + str(nEntries))

  for i in range(nEntries):
    if (i%20000 == 0):
      print(" ==> Processing " + str(i) + "th Event")
    chain.GetEntry(i)
    for sig in signals:
      #if (sig == getattr(chain,"sig")):
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

#multiSig("lxy_bs" ,50,0,3   ,"lxy_bs" ,xLabel = "[cm]")
#multiSig("lxy_ds" ,50,0,1   ,"lxy_ds" ,xLabel = "[cm]")
#multiSig("lxy_phi",50,0,1.8 ,"lxy_phi",xLabel = "[cm]")

#multiSig("dxy_mu" ,50,0,0.3 ,"dxy_mu" ,xLabel = "[cm]")
#multiSig("dxy_pi" ,50,0,0.6 ,"dxy_pi" ,xLabel = "[cm]")
#multiSig("dxy_k1" ,50,0,0.4 ,"dxy_k1" ,xLabel = "[cm]")

#multiSig("dz_mu"  ,50,0,0.3 ,"dz_mu" ,xLabel = "[cm]")
#multiSig("dz_pi"  ,50,0,0.7 ,"dz_pi" ,xLabel = "[cm]")
#multiSig("dz_k1"  ,50,0,0.5 ,"dz_k1" ,xLabel = "[cm]")

### bs momentum comparison

#multiVar(["bs_pt_coll", "bs_pt_lhcb", "bs_pt_lhcb_alt", "bs_pt_reco_1", "bs_pt_reco_2", "bs_gen_pt"], 50,0.0,150.0, "bsPtComparison")
multiVar(["q2_coll", "q2_lhcb", "q2_lhcb_alt", "q2_reco_1", "q2_reco_2", "q2_gen"], 50,0.0,20.0, "q2Comparison",xLabel = r"Q^{2} [GeV^{2}]")
multiVar(["m2_miss_coll", "m2_miss_lhcb", "m2_miss_lhcb_alt", "m2_miss_reco_1", "m2_miss_reco_2", "m2_miss_gen"], 50,0.0,10.0, "m2MissComparison", xLabel = r" m^{2}_{miss} [GeV^{2}]")

#multiVar(["bs_eta_coll", "bs_eta_lhcb", "bs_eta_lhcb_alt", "bs_eta_reco_1", "bs_eta_reco_2", "bs_gen_eta"], 50,-3,3.0, "bsEtaComparison")
#multiVar(["bs_phi_coll", "bs_phi_lhcb", "bs_phi_lhcb_alt", "bs_phi_reco_1", "bs_phi_reco_2", "bs_gen_phi"], 50,-3.2,3.2, "bsPhiComparison")
#multiVar(["bs_boost_coll", "bs_boost_lhcb", "bs_boost_lhcb_alt", "bs_boost_reco_1", "bs_boost_reco_2", "bs_boost_gen"], 50,0.9,1, "bsBoostComparison")
#multiVar(["bs_boost_coll_pt", "bs_boost_lhcb_pt", "bs_boost_lhcb_alt_pt", "bs_boost_reco_1_pt", "bs_boost_reco_2_pt", "bs_boost_gen_pt"], 50,0,1, "bsPtBoostComparison")
#multiVar(["bs_boost_coll_eta", "bs_boost_lhcb_eta", "bs_boost_lhcb_alt_eta", "bs_boost_reco_1_eta", "bs_boost_reco_2_eta", "bs_boost_gen_eta"], 50,-3,3, "bsEtaBoostComparison")
#multiVar(["bs_boost_coll_phi", "bs_boost_lhcb_phi", "bs_boost_lhcb_alt_phi", "bs_boost_reco_1_phi", "bs_boost_reco_2_phi", "bs_boost_gen_phi"], 50,-3.2,3.2, "bsPhiBoostComparison")


#hel angles

#collinear approx spoils this distribution as it ds mu system is back to back in bs rest frame --> dont plot it
#multiVar(["cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen", "cosMuWGenLhcb", "cosMuWGenReco1", "cosMuWGenReco2"], 30,-1,1, "cosMuWComparison",xLabel = "cos(#mu W)", color = colors[1:-1]) 
multiVar(["cosMuWLhcb", "cosMuWLhcbAlt", "cosMuWReco1", "cosMuWReco2", "cosMuWGen"], 30,-1,1, "cosMuWComparison",xLabel = "cos(#mu W)", color = colors[1:-1]) 

#multiVar(["cosMuWGen", "cosMuWGenLhcb", "cosMuWGenReco1", "cosMuWGenReco2"], 30,-1,1, "cosMuWGenComparison",norm = False, xLabel = "") 


multiVar(["cosPhiDsColl","cosPhiDsLhcb", "cosPhiDsLhcbAlt", "cosPhiDsReco1", "cosPhiDsReco2", "cosPhiDsGen"], 50,-1,1, "cosPhiDsComparison",xLabel = r"cos(#phi D_{s})") 
#multiVar(["cosPiDsLhcb", "cosPiDsLhcbAlt", "cosPiDsReco1", "cosPiDsReco2", "cosPiDsGen", "cosPiDsGenLhcb"], 50,-1,1, "cosPiDsComparison",xLabel = "") 
multiVar(["cosPlaneBsColl","cosPlaneBsLhcb", "cosPlaneBsLhcbAlt", "cosPlaneBsReco1", "cosPlaneBsReco2", "cosPlaneBsGen"], 30,-1,1, "cosPlaneBs",xLabel = "cos Plane 1")
multiVar(["cosPlaneDsColl","cosPlaneDsLhcb", "cosPlaneDsLhcbAlt", "cosPlaneDsReco1", "cosPlaneDsReco2", "cosPlaneDsGen"], 30,-1,1, "cosPlaneDs",xLabel = "cos Plane 2")

#multiVar(["cosPiK1","cosPiK1Gen"],30,-1,1,"cosPiK1", xLabel = "")

#compare prefit, postfit and gen momenta

#multiVar(["pv_x",    "pv_x_gen"],                     30,-0.01,0.03, "pv_x", xLabel = "[cm]", color = [ROOT.kGray + 5,ROOT.kMagenta])
#multiVar(["pv_y",    "pv_y_gen"],                     30,0.02,0.06, "pv_y", xLabel = "[cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
#multiVar(["pv_z",    "pv_z_gen"],                     30,-5,5, "pv_z", xLabel = "[cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

#multiVar(["sv_x",    "sv_x_gen"],                     30,-1,1, "sv_x", xLabel = "[cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
#multiVar(["sv_y",    "sv_y_gen"],                     30,-1,1, "sv_y", xLabel = "[cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])
#multiVar(["sv_z",    "sv_z_gen"],                     30,-5,5, "sv_z", xLabel = "[cm]",color = [ROOT.kGray + 5,ROOT.kMagenta])

#multiVar(["k1_pt",    "k1_refitted_pt",    "k1_gen_pt"],                     30,0,30, "pre_vs_postfit_k1")
#multiVar(["k2_pt",    "k2_refitted_pt",    "k2_gen_pt"],                     30,0,30, "pre_vs_postfit_k2")
#multiVar(["pi_pt",    "pi_refitted_pt",    "pi_gen_pt"],                     50,0,30, "pre_vs_postfit_pi")
#multiVar(["mu_pt",    "mu_refitted_pt",    "mu_gen_pt"],                       50,0,30, "pre_vs_postfit_mu_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])
#multiVar(["mu_eta",    "mu_refitted_eta",    "mu_gen_eta"],                    50,-3,3, "pre_vs_postfit_mu_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])
#multiVar(["mu_phi",    "mu_refitted_phi",    "mu_gen_phi"],                    50,-8,8, "pre_vs_postfit_mu_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta])

#multiVar(["kk_pt",    "phi_refitted_pt",   "phi_gen_pt", "phi_fitted_pt"],     50,0,30, "pre_vs_postfit_phi_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6 ])
#multiVar(["kk_eta",    "phi_refitted_eta",   "phi_gen_eta", "phi_fitted_eta"], 50,-3,3, "pre_vs_postfit_phi_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["kk_phi",    "phi_refitted_phi",   "phi_gen_phi", "phi_fitted_phi"], 50,-8,8, "pre_vs_postfit_phi_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])

#multiVar(["phiPi_pt", "ds_refitted_pt",    "ds_gen_pt", "ds_fitted_pt"],       50,0,30, "pre_vs_postfit_ds_pt", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["phiPi_eta", "ds_refitted_eta",    "ds_gen_eta", "ds_fitted_eta"],   50,-3,3, "pre_vs_postfit_ds_eta", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])
#multiVar(["phiPi_phi", "ds_refitted_phi",    "ds_gen_phi", "ds_fitted_phi"],   50,-8,8, "pre_vs_postfit_ds_phi", color = [ ROOT.kTeal + 6,ROOT.kAzure + 6,ROOT.kMagenta,ROOT.kViolet + 6])



