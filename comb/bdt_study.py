import ROOT

chain = ROOT.TChain("tree")
chain.Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/data_05Dec2025_11h14m48s/*0*.root")

variables = [
"mu_pt",
"k1_pt",
"k2_pt",
"pi_pt",
"bs_pt",
"kk_pt",
"phiPi_pt",

"mu_eta",
"pi_eta",
"k1_eta",
"k2_eta",

"lxy_ds_sig",
"sv_chi2",
"tv_chi2",
"fv_chi2",

"phiPi_deltaR",
"kk_deltaR",
"dsMu_deltaR",

"bs_pt_coll",
"bs_pt_lhcb_alt",
"bs_pt_reco_1",
"bs_pt_reco_2",

"dsMu_m",
"phiPi_m",

"q2_coll",
"q2_lhcb_alt",
"q2_reco_1",
"q2_reco_2",

"cosMuW_coll",
"cosMuW_lhcb_alt",
"cosMuW_reco_1",
"cosMuW_reco_2",

"e_star_coll",
"e_star_lhcb_alt",
"e_star_reco_1",
"e_star_reco_2",


]

bins = [
"(30,0,30)",
"(30,0,15)",
"(30,0,15)",
"(30,0,15)",
"(30,0,50)",
"(30,0,20)",
"(30,0,20)",

"(30,0,2.5)",
"(30,0,2.5)",
"(30,0,2.5)",
"(30,0,2.5)",

"(30,0,10)",
"(30,0,10)",
"(30,0,10)",
"(30,0,10)",

"(30,0,1.0)",
"(30,0,0.3)",
"(30,0,0.3)",

"(30,0,50)",
"(30,0,50)",
"(30,0,50)",
"(30,0,50)",

"(30,0,8)",
"(30,1.9,2.05)",

"(30,0,12)",
"(30,0,12)",
"(30,0,12)",
"(30,0,12)",

"(50,-1,1)",
"(50,-1,1)",
"(50,-1,1)",
"(50,-1,1)",

"(30,0,3)",
"(30,0,3)",
"(30,0,3)",
"(30,0,3)",

]


ROOT.gStyle.SetOptStat(0)

for (v,b) in zip(variables,bins):

  print(v)
  print(b)
  c = ROOT.TCanvas("","",800,800)
  c.SetLeftMargin(0.14)
  c.SetBottomMargin(0.12)
  c.SetTopMargin(0.06)
  c.SetRightMargin(0.04)
  chain.SetLineWidth(2)

  # kk flip inclusive
  chain.SetLineColor(ROOT.kOrange)
  chain.Draw(v + ">>h1" + b, "(k1_charge * k2_charge >0) && (dsMu_m > 5.366) ", "HIST NORM")
  h_kk = ROOT.gPad.GetPrimitive(f"h1")
  h_kk.SetTitle("")
  h_kk.GetXaxis().SetTitle(v)
  h_kk.GetYaxis().SetTitle("a.u.")


  # pimu flip
  chain.SetLineColor(ROOT.kRed)
  chain.Draw(v + ">>h2" + b, "(k1_charge * k2_charge <0) && (pi_charge*mu_charge >0) && (dsMu_m > 5.366) ", "HIST NORM, SAME")
  h_pimu = ROOT.gPad.GetPrimitive(f"h2")

  #correct sign events
  chain.SetLineColor(ROOT.kGreen)
  chain.Draw(v + ">>h3" + b, "((pi_charge*mu_charge < 0) && (k1_charge * k2_charge <0)) && (dsMu_m > 5.366)", "HIST NORM SAME")
  h_corr = ROOT.gPad.GetPrimitive(f"h3")

  max_y = max(h_kk.GetMaximum(), h_pimu.GetMaximum(), h_corr.GetMaximum())
  h_kk.SetMaximum(1.3 * max_y)

  h_kk.Draw("HIST")
  h_pimu.Draw("HIST SAME")
  h_corr.Draw("HIST SAME")

  # legend
  leg = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)

  leg.AddEntry(h_kk,   "KK wrong", "l")
  leg.AddEntry(h_pimu, "#pi #mu wrong", "l")
  leg.AddEntry(h_corr, "Correct Sign", "l")
  leg.Draw()

  c.SaveAs(f"./bdt_study/{v}.pdf" )
  
