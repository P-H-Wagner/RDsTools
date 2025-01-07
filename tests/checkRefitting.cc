#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}


void print(vector<string> vars, ROOT::RDataFrame d, vector<EColor> col, string signal, string name){

  TCanvas* c = new TCanvas("", "Canvas", 800, 600);
  c->cd();
  TLegend leg = getLegend();
  
  const char* v0 = vars.at(0).c_str();
  const char* v1 = vars.at(1).c_str();
  const char* v2 = vars.at(2).c_str();

  int bins; double start; double stop;

  bins = 20;

  if (name.find("pt") != string::npos){
    start = 7;
    stop = 12;
  }
  else if (name.find("eta") != string::npos){
    start = -2.5;
    stop = 2.5;
  }
  else{
    start = -3.14;
    stop = +3.14;
  }


  auto h0 = d.Filter(signal).Histo1D({"","",bins,start,stop}, v0);
  auto h1 = d.Filter(signal).Histo1D({"","",bins,start,stop}, v1);
  auto h2 = d.Filter(signal).Histo1D({"","",bins,start,stop}, v2);

  h0->SetLineColor(col.at(0));
  h1->SetLineColor(col.at(1));
  h2->SetLineColor(col.at(2));

  leg.AddEntry(h0.GetPtr(), v0, "L");
  leg.AddEntry(h1.GetPtr(), v1, "L");
  leg.AddEntry(h2.GetPtr(), v2, "L");

  h0->Draw("HIST");
  h1->Draw("HIST SAME");
  h2->Draw("HIST SAME");

  string toSave = name + ".pdf";
  cout << "saving to... " << toSave << endl;
  c->SaveAs(toSave.c_str());

}


void checkRefitting(){

  // take one signal file as example
  const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/26_07_2024_14_46_03/all_signals_flatChunk_1*";
 
  //select signal
  int sig = 1;
  map<int, string> folder;

  folder[0]  = "dsmu";
  folder[1]  = "dstau";
  folder[10] = "dsstarmu";
  folder[11] = "dsstartau";

  string signal = "(gen_sig == " + to_string(sig) + ")";

  //load trees
  ROOT::RDataFrame df("tree", fin);


  vector<string> vars_mu_pt;
  vector<string> vars_mu_eta;
  vector<string> vars_mu_phi;

  vector<string> vars_pi_pt;
  vector<string> vars_pi_eta;
  vector<string> vars_pi_phi;

  vector<string> vars_k1_pt;
  vector<string> vars_k1_eta;
  vector<string> vars_k1_phi;

  vector<string> vars_k2_pt;
  vector<string> vars_k2_eta;
  vector<string> vars_k2_phi;

  vars_mu_pt.push_back("gen_mu_pt");
  vars_mu_pt.push_back("mu_pt");
  vars_mu_pt.push_back("mu_refitted_pt");

  vars_mu_eta.push_back("gen_mu_eta");
  vars_mu_eta.push_back("mu_eta");
  vars_mu_eta.push_back("mu_refitted_eta");

  vars_mu_phi.push_back("gen_mu_phi");
  vars_mu_phi.push_back("mu_phi");
  vars_mu_phi.push_back("mu_refitted_phi");

  vars_pi_pt.push_back("gen_pi_pt");
  vars_pi_pt.push_back("pi_pt");
  vars_pi_pt.push_back("pi_refitted_pt");

  vars_pi_eta.push_back("gen_pi_eta");
  vars_pi_eta.push_back("pi_eta");
  vars_pi_eta.push_back("pi_refitted_eta");

  vars_pi_phi.push_back("gen_pi_phi");
  vars_pi_phi.push_back("pi_phi");
  vars_pi_phi.push_back("pi_refitted_phi");

  vars_k1_pt.push_back("gen_k1_pt");
  vars_k1_pt.push_back("k1_pt");
  vars_k1_pt.push_back("k1_refitted_pt");
  
  vars_k1_eta.push_back("gen_k1_eta");
  vars_k1_eta.push_back("k1_eta");
  vars_k1_eta.push_back("k1_refitted_eta");

  vars_k1_phi.push_back("gen_k1_phi");
  vars_k1_phi.push_back("k1_phi");
  vars_k1_phi.push_back("k1_refitted_phi");

  vars_k2_pt.push_back("gen_k2_pt");
  vars_k2_pt.push_back("k2_pt");
  vars_k2_pt.push_back("k2_refitted_pt");

  vars_k2_eta.push_back("gen_k2_eta");
  vars_k2_eta.push_back("k2_eta");
  vars_k2_eta.push_back("k2_refitted_eta");

  vars_k2_phi.push_back("gen_k2_phi");
  vars_k2_phi.push_back("k2_phi");
  vars_k2_phi.push_back("k2_refitted_phi");

  vector<EColor> colors;
  colors.push_back(kPink);
  colors.push_back(kBlack);
  colors.push_back(kBlue);
 
  print(vars_mu_pt,  df, colors, signal, "./" + folder[sig] + "/mu_pt");
  print(vars_pi_pt,  df, colors, signal, "./" + folder[sig] + "/pi_pt");
  print(vars_k1_pt,  df, colors, signal, "./" + folder[sig] + "/k1_pt");
  print(vars_k2_pt,  df, colors, signal, "./" + folder[sig] + "/k2_pt");

  print(vars_mu_eta, df, colors, signal, "./" + folder[sig] + "/mu_eta");
  print(vars_pi_eta, df, colors, signal, "./" + folder[sig] + "/pi_eta");
  print(vars_k1_eta, df, colors, signal, "./" + folder[sig] + "/k1_eta");
  print(vars_k2_eta, df, colors, signal, "./" + folder[sig] + "/k2_eta");

  print(vars_mu_phi, df, colors, signal, "./" + folder[sig] + "/mu_phi");
  print(vars_pi_phi, df, colors, signal, "./" + folder[sig] + "/pi_phi");
  print(vars_k1_phi, df, colors, signal, "./" + folder[sig] + "/k1_phi");
  print(vars_k2_phi, df, colors, signal, "./" + folder[sig] + "/k2_phi");

}
