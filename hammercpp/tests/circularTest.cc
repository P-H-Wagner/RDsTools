#include "Hammer/Tools/HammerRoot.hh"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h>

using namespace std;

// define input files
const char* getInputFile(string signal){

  const char* fin;

  if (signal == "dsmu") fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/gen/dsmu_ISGW2_20_05_2025_14_37_46/*"; // old: /dsmu_ISGW2_10_12_2024_20_00_43/*"; //produced with HQET2
  else fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/gen/dsmu_isgw2_CLN_20_05_2025_14_32_36/*"; // old: dsmu_isgw2_CLN_10_12_2024_19_31_12/*"; //produced with ISGW2

  return fin;

}

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}
int main(){

  gStyle->SetOptStat(0);

  //load hammer trees
  ROOT::RDataFrame df("tree",  getInputFile("dsmu") );  
  ROOT::RDataFrame df_isgw2("tree",  getInputFile("dsmu_isgw2") );  

  //first fill dsmu
  auto h         = df.Filter("(gen_sig == 0)").Histo1D({"h", "", 20, 0, 12}, "gen_q2");
  auto h_w       = df.Filter("(gen_sig == 0)").Histo1D({"h", "", 20, 0, 12}, "gen_q2", "central_w");

  //normalize
  h->Scale(1/h->Integral());
  h_w->Scale(1/h_w->Integral());

  h->SetLineColor(kBlue);
  h_w->SetLineColor(kBlue);
  h_w->SetLineStyle(7);
  
  //now fill isgw2
  auto h_isgw2   = df_isgw2.Filter("(gen_sig == 0)").Histo1D({"h_isgw2", "", 20, 0, 12}, "gen_q2");
  auto h_isgw2_w = df_isgw2.Filter("(gen_sig == 0)").Histo1D({"h_isgw2", "", 20, 0, 12}, "gen_q2", "central_w");

  //cout << "#event simulated with heqt2: " << h->GetEntries() << endl;
  //cout << "#event simulated with isgw2: " << h_isgw2->GetEntries() << endl;

  h_isgw2->Scale(1/h_isgw2->Integral());
  h_isgw2_w->Scale(1/h_isgw2_w->Integral());


  h_isgw2->SetLineColor(kRed);
  h_isgw2_w->SetLineColor(kRed);
  h_isgw2_w->SetLineStyle(7);

  ////////////////////////
  // Draw unweighted    //
  ////////////////////////

  TCanvas *c1 = new TCanvas("c1", "My Canvas", 800, 600);
  
  h->SetMaximum(0.12); // Set y-axis maximum
  h->GetXaxis()->SetTitle("gen. Q^{2} (GeV^{2})");
  h->GetYaxis()->SetTitle("a.u.");

  h->Draw("HIST ");
  h_isgw2->Draw("HIST SAME");

  TLegend leg = getLegend();
  leg.AddEntry(h.GetPtr(), "HQET2", "L");
  leg.AddEntry(h_isgw2.GetPtr(), "ISGW2", "L");

  leg.Draw("SAME");
  

  c1->Update();
  c1->SaveAs("wout_weights.pdf");

  //////////////////////////////////
  // Draw hqet2 vs. weighted hqet2//
  //////////////////////////////////

  TCanvas *c2 = new TCanvas("c1", "My Canvas", 800, 600);
  
  h->SetMaximum(0.12); // Set y-axis maximum
  h->GetXaxis()->SetTitle("gen. Q^{2} (GeV^{2})");
  h->GetYaxis()->SetTitle("a.u.");

  h->Draw("HIST ");
  h_w->Draw("HIST SAME");

  TLegend leg2 = getLegend();
  leg2.AddEntry(h.GetPtr(), "HQET2", "L");
  leg2.AddEntry(h_w.GetPtr(), "Weighted HQET2", "L");

  leg2.Draw("SAME");
  

  c2->Update();
  c2->SaveAs("weighted_hqet2.pdf");

  /////////////////////////////
  // Now weigh hqet2 -> isgw2//
  /////////////////////////////

  TCanvas *c3 = new TCanvas("c3", "My Canvas", 800, 600);
  
  h_isgw2->SetMaximum(0.12); // Set y-axis maximum
  h_isgw2->GetXaxis()->SetTitle("gen. Q^{2} (GeV^{2})");
  h_isgw2->GetYaxis()->SetTitle("a.u.");


  h_isgw2->Draw("HIST ");
  h_w->Draw("HIST SAME");

  TLegend leg3 = getLegend();
  leg3.AddEntry(h_w.GetPtr(), "Weighted HQET2", "L");
  leg3.AddEntry(h_isgw2.GetPtr(), "ISGW2", "L");

  leg3.Draw("SAME");
  

  c3->Update();
  c3->SaveAs("hqet2_to_isgw2.pdf");

  /////////////////////////////
  // Now weigh isgw2 -> hqet2//
  /////////////////////////////

  TCanvas *c4 = new TCanvas("c4", "My Canvas", 800, 600);
  
  h->SetMaximum(0.12); // Set y-axis maximum
  h->GetXaxis()->SetTitle("gen. Q^{2} (GeV^{2})");
  h->GetYaxis()->SetTitle("a.u.");
  h->Draw("HIST ");
  h_isgw2_w->Draw("HIST SAME");

  TLegend leg4 = getLegend();
  leg4.AddEntry(h.GetPtr(), "HQET2", "L");
  leg4.AddEntry(h_isgw2_w.GetPtr(), "Weighted ISGW2", "L");

  leg4.Draw("SAME");
  

  c4->Update();
  c4->SaveAs("isgw2_to_hqet2.pdf");







  /*
  // now weigh isgw2 -> CLN 
  TCanvas *c2 = new TCanvas("c1", "My Canvas", 800, 600);
  tree->SetLineColor(kBlue);
  tree->Draw("gen_q2>>h",sel ,"HIST NORM");
  tree_isgw2->SetLineColor(kRed);
  tree_isgw2->Draw("gen_q2", "central_w * (" + sel + ")"  ,"HIST NORM SAME");

  h->SetMinimum(0);  // Set y-axis minimum
  h->SetMaximum(0.12); // Set y-axis maximum


  c2->SaveAs("weigh_isgw2.pdf");

  // now weigh CLN -> isgw2 
  TCanvas *c3 = new TCanvas("c1", "My Canvas", 800, 600);

  tree_isgw2->SetLineColor(kRed);
  tree_isgw2->Draw("gen_q2>>h", sel ,"HIST NORM");
  tree->SetLineColor(kBlue);
  tree->Draw("gen_q2",  "central_w * (" + sel + ")" ,"HIST NORM SAME");

  h->SetMinimum(0);  // Set y-axis minimum
  h->SetMaximum(0.12); // Set y-axis maximum


  c3->SaveAs("weigh_cln.pdf");

  TCanvas *c4 = new TCanvas("c1", "My Canvas", 800, 600);

  tree_isgw2->SetLineColor(kRed);
  tree_isgw2->Draw("gen_q2>>h", sel ,"HIST NORM");
  tree->SetLineColor(kBlue);
  tree->Draw("gen_q2",  "central_w * (" + sel + ")" ,"HIST NORM SAME");

  h->SetMinimum(0);  // Set y-axis minimum
  h->SetMaximum(0.12); // Set y-axis maximum


  c3->SaveAs("weigh_cln.pdf");
  */

}


