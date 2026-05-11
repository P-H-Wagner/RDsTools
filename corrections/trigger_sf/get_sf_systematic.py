import ROOT
import numpy as np
import argparse
from datetime import datetime
import os
import sys
sys.path.append(os.path.abspath("/work/pahwagne/RDsTools/help"))
from helper import *

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)


selection = minimal

#variables for which we want to check the shape effect
var = ["q2_coll","phiPi_m","score1"]

models = {}
models["q2_coll"] = ROOT.RDF.TH1DModel("q2_coll", '',   31,       0,     12)
models["phiPi_m"] = ROOT.RDF.TH1DModel("phiPi_m", '',   31,    1.91,  2.028) 
models["score1"]  = ROOT.RDF.TH1DModel("score1" , '',   31,       0,    0.8)

#import sf and edges
scale_factors     = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_fullBPark.txt")
scale_factors_err = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/trigger_sf_err_fullBPark.txt")
pt_edges          = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/pt_edges_fullBPark.txt")
dxy_edges         = np.loadtxt("/work/pahwagne/RDsTools/corrections/trigger_sf/dxy_edges_fullBPark.txt")

#convert into string such that it can be handed to the cpp process line
pt_edges_str  = ", ".join(map(str, pt_edges     )) # of type: '0.01, 0.34, 0.89'
dxy_edges_str = ", ".join(map(str, dxy_edges    )) # of type: '0.01, 0.34, 0.89'

#scale_factors is a matrix!
lenx = scale_factors.shape[0]
leny = scale_factors.shape[1]

#create RDF 
#input file is a gen lvl sample wout filter

chain_sig = ROOT.TChain("tree")
#chain_sig.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_{nn_model}/*1*")
chain_sig.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano//skimmed/minimal/17_03_2026_11_31_56/*1*")

df_sig = ROOT.RDataFrame(chain_sig)

chain_hb = ROOT.TChain("tree")
#chain_hb.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/hb_{nn_model}/*1*")
chain_hb.Add(f"/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano//skimmed/minimal/25_02_2026_16_36_40/*1*")
df_hb = ROOT.RDataFrame(chain_hb)

now = datetime.now()
dt  = now.strftime("%d_%m_%Y_%H_%M_%S")

histos_sig     = {}
histos_sig_err_up = {}
histos_sig_err_down = {}
histos_hb      = {}
histos_hb_err_up  = {}
histos_hb_err_down  = {}


#we loop over all the pt,eta bins

for v in var:

  #where to save the plots for var v 
  dest_dir = f"/work/pahwagne/RDsTools/corrections/trigger_sf/nn_model_{nn_model}/{v}/"
  os.system(f"mkdir -p {dest_dir}")


  #now loop over all the sf bins
  for i in range(len(pt_edges)-1): #array([  6. ,   7. ,   8. ,   8.5,   9. ,  10. ,  10.5,  11. ,  12. , 20. , 100. ])

    print(f"====> Pt bin {pt_edges[i]} < pt < {pt_edges[i+1]}")

    for j in range(len(dxy_edges)-1): #array([  0. ,   3. ,   3.5,   4. ,   5. ,   6. ,   8. ,  10. ,  20. , 500. ]) 

      print(f"====> dxy bin {dxy_edges[j]} < dxy < {dxy_edges[j+1]}")
 
      print(f"Scale factor is: {scale_factors    [i][j]}")
      print(f"err is         : {scale_factors_err[i][j]}")

                            #bool to check if we are in this pt bin
      sel_str          = f" ( abs(trigger_sf - {scale_factors[i][j]}) < 1e-6 ) * trigger_sf"         
      sel_str_err_up   = f"(( abs(trigger_sf - {scale_factors[i][j]}) < 1e-6 ) * (trigger_sf + trigger_sf_err)) + (( abs(trigger_sf - {scale_factors[i][j]}) > 1e-6 ) * (trigger_sf )) "         
      sel_str_err_down = f"(( abs(trigger_sf - {scale_factors[i][j]}) < 1e-6 ) * (trigger_sf - trigger_sf_err)) + (( abs(trigger_sf - {scale_factors[i][j]}) > 1e-6 ) * (trigger_sf )) "         
    
      print(sel_str)
      print(sel_str_err_up)

      #Fill and plot histo!!

      #fill the histogram vor var v

      #apply the central trigger scale factors (always)
      histos_sig[v]     = df_sig.Filter(selection).Histo1D(models[v], v , "trigger_sf")
      histos_sig[v].SetLineColor(ROOT.kBlue)    
      histos_hb[v]      = df_hb .Filter(selection).Histo1D(models[v], v , "trigger_sf")
      histos_hb[v].SetLineColor(ROOT.kBlue)    
    
      #w.l.o.g. we can only consider the up variation (down err is the same!)

      #apply trigger sf variation only in bin ij, otherwise the central one
      histos_sig_err_up[v] = df_sig.Filter(selection).Define("bin_ij", sel_str_err_up).Histo1D(models[v], v , "bin_ij")
      histos_sig_err_up[v].SetLineColor(ROOT.kRed)    
      histos_sig_err_up[v].SetLineStyle(7)    

      histos_hb_err_up[v]  = df_hb .Filter(selection).Define("bin_ij", sel_str_err_up).Histo1D(models[v], v , "bin_ij")
      histos_hb_err_up[v].SetLineColor(ROOT.kRed)    
      histos_hb_err_up[v].SetLineStyle(7)    

      histos_sig_err_down[v] = df_sig.Filter(selection).Define("bin_ij", sel_str_err_down).Histo1D(models[v], v , "bin_ij")
      histos_sig_err_down[v].SetLineColor(ROOT.kGreen)    
      histos_sig_err_down[v].SetLineStyle(4)    

      histos_hb_err_down[v]  = df_hb .Filter(selection).Define("bin_ij", sel_str_err_down).Histo1D(models[v], v , "bin_ij")
      histos_hb_err_down[v].SetLineColor(ROOT.kGreen)    
      histos_hb_err_down[v].SetLineStyle(4)    

      canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
      canvas.Draw()
      canvas.cd()

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)


      histos_sig[v].GetYaxis().SetRangeUser(1e-3, histos_sig[v].GetMaximum() * 1.5 )
      histos_sig[v].Draw("HIST")
      histos_sig_err_up[v].Draw("HIST SAME")
      histos_sig_err_down[v].Draw("HIST SAME")

      legend.AddEntry(histos_sig[v].GetPtr()    ,"central"   ,"L")
      legend.AddEntry(histos_sig_err_up[v]  .GetPtr(),r"#sigma up","L")
      legend.AddEntry(histos_sig_err_down[v].GetPtr(),r"#sigma down","L")
      legend.Draw("SAME")     
 
      canvas.SaveAs(f"{dest_dir}/bin_pt_bin_{i}_dxy_bin_{j}_sig.pdf")

      canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
      canvas.Draw()
      canvas.cd()

      legend = ROOT.TLegend(0.5,0.8,0.9,0.9)
      legend.SetTextSize(0.04)
      legend.SetBorderSize(0)
      legend.SetFillStyle(0)


      histos_hb[v].GetYaxis().SetRangeUser(1e-3, histos_hb[v].GetMaximum() * 1.5 )
      histos_hb[v].Draw("HIST")
      histos_hb_err_up[v].Draw("HIST SAME")
      histos_hb_err_down[v].Draw("HIST SAME")

      legend.AddEntry(histos_hb[v].GetPtr()    ,"central"   ,"L")
      legend.AddEntry(histos_hb_err_up[v]  .GetPtr(),r"#hbma up","L")
      legend.AddEntry(histos_hb_err_down[v].GetPtr(),r"#hbma down","L")
      legend.Draw("SAME")     
 
      canvas.SaveAs(f"{dest_dir}/bin_pt_bin_{i}_dxy_bin_{j}_hb.pdf")


