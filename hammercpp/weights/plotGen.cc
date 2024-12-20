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

  if (signal == "dsstarmu") {
    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/dsstarmu_BGLVar_18_12_2024_10_36_23/*"; //first harrison production
    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/dsstarmu_BGLVar_20_12_2024_14_25_22/*"; //second harrison with F1_0 set
  }
  else {
    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/dsstartau_BGLVar_18_12_2024_10_36_07/*"; //produced with ISGW2
  }
  return fin;

}

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}

void plot(ROOT::RDataFrame df, int bins,double start,double stop, string var){


  ROOT::RDF::RResultPtr<TH1D> h;
  ROOT::RDF::RResultPtr<TH1D> h_w;

  if (var.find("eta") !=  std::string::npos) {

    h   = df.Filter("(gen_sig == 10)").Define("abs_eta", [](double x) { return std::abs(x); }, {var}).Histo1D({"h", "", bins, start, stop}, "abs_eta");
    h_w = df.Filter("(gen_sig == 10)").Define("abs_eta", [](double x) { return std::abs(x); }, {var}).Histo1D({"h", "", bins, start, stop}, "abs_eta", "central_w");

  }
  else{
  
    h   = df.Filter("(gen_sig == 10)").Histo1D({"h", "", bins, start, stop}, var);
    h_w = df.Filter("(gen_sig == 10)").Histo1D({"h", "", bins, start, stop}, var, "central_w");
  }

  h->SetLineColor(kRed);
  h_w->SetLineColor(kBlue);

  h->Scale(1/h->Integral());
  h_w->Scale(1/h_w->Integral());

  if (var.find("pt") != std::string::npos){
  
    h->GetXaxis()->SetRangeUser(7.0,12);
    h_w->GetXaxis()->SetRangeUser(7.0,12);

  }

  TCanvas *c1 = new TCanvas("c1", "My Canvas", 800, 600);
  
  h->SetMaximum(h->GetMaximum() * 1.3); // Set y-axis maximum
  h->Draw("HIST ");
  h_w->Draw("HIST SAME");

  TLegend leg = getLegend();
  leg.AddEntry(h.GetPtr(), "Unweighted", "L");
  leg.AddEntry(h_w.GetPtr(), "Weighted", "L");

  leg.Draw("SAME");
  

  c1->Update();
  string toSave = "./gen/" + var + ".pdf";
  c1->SaveAs(toSave.c_str());

  for(size_t i = 1; i < 11; ++i){

    auto h_up   = df.Filter("(gen_sig == 10)").Histo1D({"h", "", bins, start, stop}, var, "e"+to_string(i)+"_up");
    auto h_down = df.Filter("(gen_sig == 10)").Histo1D({"h", "", bins, start, stop}, var, "e"+to_string(i)+"_down");

    h_up->SetLineColor(kBlack);
    h_down->SetLineColor(kBlack);

    h_up->Scale(1.0/h_up->Integral());
    h_down->Scale(1.0/h_down->Integral());

    if (var.find("pt") != std::string::npos){
    
      h_up->GetXaxis()->SetRangeUser(7.0,12);
      h_down->GetXaxis()->SetRangeUser(7.0,12);
  
    }

    TCanvas *c2 = new TCanvas("c2", "My Canvas", 800, 600);
    
    h_w->SetMaximum(h->GetMaximum() * 1.3); // Set y-axis maximum
    h_w->Draw("HIST ");
    h_up->Draw("HIST SAME");
    h_down->Draw("HIST SAME");
  
    TLegend leg2 = getLegend();
    leg2.AddEntry(h_w.GetPtr(), "Central", "L");
    leg2.AddEntry(h_up.GetPtr(), "Up", "L");
    leg2.AddEntry(h_down.GetPtr(), "Down", "L");
  
    leg2.Draw("SAME");
    
  
    c2->Update();
    string toSave = "./gen/" + var + "_var_" + to_string(i) + ".pdf";
    c2->SaveAs(toSave.c_str());
  





  }

}

int main(){

  gStyle->SetOptStat(0);

  //load hammer trees
  ROOT::RDataFrame df("tree",  getInputFile("dsstarmu") );  
  //ROOT::RDataFrame df_isgw2("tree",  getInputFile("dsstartau") );  

  //first fill dsmu
  plot(df, 20, 0, 12, "gen_mu_pt");
  plot(df, 20, 0, 12, "gen_pi_pt");
  plot(df, 20, 0, 12, "gen_k1_pt");
  plot(df, 20, 0, 12, "gen_k2_pt");
  plot(df, 10, 0, 5, "gen_mu_eta");
  plot(df, 20, 0, 5, "gen_pi_eta");
  plot(df, 20, 0, 5, "gen_k1_eta");
  plot(df, 20, 0, 5, "gen_k2_eta");
  plot(df, 20, -3.14, 3.14, "gen_mu_phi");
  plot(df, 20, -3.14, 3.14, "gen_pi_phi");
  plot(df, 20, -3.14, 3.14, "gen_k1_phi");
  plot(df, 20, -3.14, 3.14, "gen_k2_phi");
  //auto h_w       = df.Filter("(gen_sig == 0)").Histo1D({"h", "", 20, 0, 12}, "gen_q2", "central_w");


}


