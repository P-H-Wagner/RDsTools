#include "Hammer/Tools/HammerRoot.hh"
#include "TTree.h"
#include "TLatex.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h>
#include <yaml-cpp/yaml.h>
#include <fstream> 

////////////////////////////////////////////////////////////////////////////////////////////
// Command line args:                                                                     //
// args[1] = BGLVar, ... (i.e. the target we want to plot with weights)                   //
// args[2] = paper (i.e. the paper we used, f.e. cohen, harrison, ... )                   //
////////////////////////////////////////////////////////////////////////////////////////////


using namespace std;

// define input files
string getInputFile(string signal, string paper){

  string fin_str;
  if (paper == "cohen"){
  fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_BGLVar_13_12_2024_13_28_24/*";
  }
  else {
  fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_BGLVar_18_12_2024_10_35_08/*";
  }
  return fin_str;

}

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}


void plotWeightEffect(auto h_wout, auto h, string key){

  //get yield effect
  double n = h->Integral();
  double n_wout = h_wout->Integral();
  int diff = int(abs(n-n_wout)/n * 100); //*100 to get %
  string yield = to_string(diff) + "% change in yield";

  const char* decay = "";
  if (key.find("dsmu")      != string::npos ){decay = "B_{s} #rightarrow D_{s} #mu #nu";}
  if (key.find("dstau")     != string::npos ){decay = "B_{s} #rightarrow D_{s} #tau #nu";}
  if (key.find("dsstarmu")  != string::npos ){decay = "B_{s} #rightarrow D*_{s} #mu #nu";}
  if (key.find("dsstartau") != string::npos ){decay = "B_{s} #rightarrow D*_{s} #tau #nu";}
  TText* text         = new TLatex(0.2, 0.85, decay);
  text->SetNDC();   
  text->SetTextFont(43);
  text->SetTextSize(25);

  TText* yield_change = new TLatex(0.5, 0.75, yield.c_str());
  yield_change->SetNDC();   
  yield_change->SetTextFont(43);
  yield_change->SetTextSize(25);

  TCanvas *c = new TCanvas("canvas", "Canvas Title", 800, 600);

  h->SetLineStyle(7);
 
  h_wout->SetMaximum(h_wout->GetMaximum()*1.4); 
  h_wout->GetYaxis()->SetTitle("events");
  h_wout->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");

  h_wout->Draw("HIST");
  h->Draw("HIST SAME");

  TLegend leg = getLegend();
  leg.AddEntry(h_wout.GetPtr(), "Unweighted", "L");
  leg.AddEntry(h.GetPtr(), "Weighted","L");
  leg.Draw("SAME");
  text->Draw("SAME");
  yield_change->Draw("SAME");

  const char* toSave = key.c_str();
  c->Update();
  c->SaveAs(toSave);


}


int main(int nargs, char* args[]){

  gStyle->SetOptStat(0);

  //variable for which we want to compute the average weight
  //string var(args[2]); //f.e. q2_coll
  string var = "class";

  //load average weights and settings 

  YAML::Node average_weights;

  if (strcmp(args[2], "cohen") == 0) {
  average_weights = YAML::LoadFile("average_weights.yaml"); //from getAverageWeights.cc
  }
  else{
  average_weights = YAML::LoadFile("average_weightsharrison.yaml"); //from getAverageWeights.cc
  }

  YAML::Node settings        = YAML::LoadFile("average_models.yaml"); //from histModels.py
  int bins   = settings[var]["bins"].as<int>();
  double min = settings[var]["xmin"].as<double>();
  double max = settings[var]["xmax"].as<double>();

  ROOT::RDataFrame* df = nullptr;

  try{
    df = new ROOT::RDataFrame("tree",  getInputFile(args[1] , string(args[2])  ));  
  }
  catch(const exception& e){ cout << "no file found" << endl; exit(1); }


  auto df_dsMu       = df->Filter("gen_sig == 0");
  auto df_dsTau      = df->Filter("gen_sig == 1");
  auto df_dsStarMu   = df->Filter("gen_sig == 10");
  auto df_dsStarTau  = df->Filter("gen_sig == 11");

  auto h_dsMu_wout        = df_dsMu.Histo1D(     {"h_dsMu",      "", bins, min, max}, var);
  auto h_dsTau_wout       = df_dsTau.Histo1D(    {"h_dsTau",     "", bins, min, max}, var);
  auto h_dsStarMu_wout    = df_dsStarMu.Histo1D( {"h_dsStarMu",  "", bins, min, max}, var);
  auto h_dsStarTau_wout   = df_dsStarTau.Histo1D({"h_dsStarTau", "", bins, min, max}, var);

  auto h_dsMu        = df_dsMu.Histo1D(     {"h_dsMu",      "", bins, min, max}, var, "central_w");
  auto h_dsTau       = df_dsTau.Histo1D(    {"h_dsTau",     "", bins, min, max}, var, "central_w");
  auto h_dsStarMu    = df_dsStarMu.Histo1D( {"h_dsStarMu",  "", bins, min, max}, var, "central_w");
  auto h_dsStarTau   = df_dsStarTau.Histo1D({"h_dsStarTau", "", bins, min, max}, var, "central_w");

  cout << "average weight is:" <<  average_weights["central_w_dsstarmu"].as<double>() << endl;
  h_dsMu->Scale(     1 / average_weights["central_w_dsmu"].as<double>());
  h_dsTau->Scale(    1 / average_weights["central_w_dstau"].as<double>());
  h_dsStarMu->Scale( 1 / average_weights["central_w_dsstarmu"].as<double>());
  h_dsStarTau->Scale(1 / average_weights["central_w_dsstartau"].as<double>());


  //First, see the effect of the central weight
  plotWeightEffect(h_dsMu_wout,      h_dsMu ,      "dsmu_weight_effect" + string(args[2]) + ".pdf");
  plotWeightEffect(h_dsTau_wout,     h_dsTau ,     "dstau_weight_effect" + string(args[2]) + ".pdf");
  plotWeightEffect(h_dsStarMu_wout,  h_dsStarMu ,  "dsstarmu_weight_effect" + string(args[2]) + ".pdf");
  plotWeightEffect(h_dsStarTau_wout, h_dsStarTau , "dsstartau_weight_effect" + string(args[2]) + ".pdf");


  h_dsMu->SetLineColor(kBlack); 
  h_dsTau->SetLineColor(kBlack); 
  h_dsStarMu->SetLineColor(kBlack); 
  h_dsStarTau->SetLineColor(kBlack); 

  vector<string> directions = {"up","down"};

  for (size_t i = 1; i < 11; i++){

    string key_dsMu      = "e"+to_string(i)+"_dsmu";
    string key_dsTau     = "e"+to_string(i)+"_dstau";
    string key_dsStarMu  = "e"+to_string(i)+"_dsstarmu";
    string key_dsStarTau = "e"+to_string(i)+"_dsstartau";

    string toSave_dsMu      = key_dsMu + string(args[2]) + ".pdf";
    string toSave_dsTau     = key_dsTau + string(args[2]) + ".pdf";
    string toSave_dsStarMu  = key_dsStarMu + string(args[2]) + ".pdf";
    string toSave_dsStarTau = key_dsStarTau + string(args[2]) + ".pdf";



    auto h_up_dsMu       = df_dsMu.Histo1D(     {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsTau      = df_dsTau.Histo1D(    {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsStarMu   = df_dsStarMu.Histo1D( {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsStarTau  = df_dsStarTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 

    h_up_dsMu->Scale(      1.0 / average_weights[ "e"+to_string(i)+"_up_dsmu"].as<double>()); 
    h_up_dsTau->Scale(     1.0 / average_weights[ "e"+to_string(i)+"_up_dstau"].as<double>()); 
    h_up_dsStarMu->Scale(  1.0 / average_weights[ "e"+to_string(i)+"_up_dsstarmu"].as<double>()); 
    h_up_dsStarTau->Scale( 1.0 / average_weights[ "e"+to_string(i)+"_up_dsstartau"].as<double>()); 

    h_up_dsMu->SetLineColor(kRed); 
    h_up_dsTau->SetLineColor(kRed); 
    h_up_dsStarMu->SetLineColor(kRed); 
    h_up_dsStarTau->SetLineColor(kRed); 

    auto h_down_dsMu       = df_dsMu.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsTau      = df_dsTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsStarMu   = df_dsStarMu.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsStarTau  = df_dsStarTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 

    h_down_dsMu->Scale(     1.0 / average_weights[ "e"+to_string(i)+"_down_dsmu"].as<double>()); 
    h_down_dsTau->Scale(    1.0 / average_weights[ "e"+to_string(i)+"_down_dstau"].as<double>()); 
    h_down_dsStarMu->Scale( 1.0 / average_weights[ "e"+to_string(i)+"_down_dsstarmu"].as<double>()); 
    h_down_dsStarTau->Scale(1.0 / average_weights[ "e"+to_string(i)+"_down_dsstartau"].as<double>()); 

    h_down_dsMu->SetLineColor(kRed); 
    h_down_dsTau->SetLineColor(kRed); 
    h_down_dsStarMu->SetLineColor(kRed); 
    h_down_dsStarTau->SetLineColor(kRed); 




    h_up_dsMu->SetMaximum(1.3 * h_up_dsMu->GetMaximum()); 
    h_up_dsTau->SetMaximum(1.3 * h_up_dsTau->GetMaximum()); 
    h_up_dsStarMu->SetMaximum(1.3 * h_up_dsStarMu->GetMaximum()); 
    h_up_dsStarTau->SetMaximum(1.3 * h_up_dsStarTau->GetMaximum()); 

    TCanvas *canvMu = new TCanvas("canvas", "Canvas Title", 800, 600);

    h_up_dsMu->Draw("HIST");
    h_dsMu->Draw("HIST SAME");
    h_down_dsMu->Draw("HIST SAME");
    canvMu->Update();
    canvMu->SaveAs(toSave_dsMu.c_str());

    cout << " integral  " <<  h_up_dsMu->Integral();
    cout << " integral  " <<  h_dsMu->Integral();
    cout << " integral  " <<  h_down_dsMu->Integral();

    TCanvas *canvTau = new TCanvas("canvas", "Canvas Title", 800, 600);
    h_up_dsTau->Draw("HIST");
    h_dsTau->Draw("HIST SAME");
    h_down_dsTau->Draw("HIST SAME");
    canvTau->Update();
    canvTau->SaveAs(toSave_dsTau.c_str());

    TCanvas *canvMuStar = new TCanvas("canvas", "Canvas Title", 800, 600);

    TLegend legMuStar = getLegend();
    legMuStar.AddEntry(h_up_dsStarMu.GetPtr(), "#pm 1#sigma envelope", "L");
    legMuStar.AddEntry(h_dsStarMu.GetPtr(), "central value", "L");

    h_up_dsStarMu->GetYaxis()->SetTitle("events");
    h_up_dsStarMu->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");

    TText* textMuStar         = new TLatex(0.2, 0.85, "B_{s} #rightarrow D*_{s} #mu #nu");
    textMuStar->SetNDC();   
    textMuStar->SetTextFont(43);
    textMuStar->SetTextSize(25);

    h_up_dsStarMu->Draw("HIST");
    h_dsStarMu->Draw("HIST SAME");
    h_down_dsStarMu->Draw("HIST SAME");
    legMuStar.Draw("SAME");
    textMuStar->Draw("SAME");

    canvMuStar->Update();
    canvMuStar->SaveAs(toSave_dsStarMu.c_str());


    TCanvas *canvTauStar = new TCanvas("canvas", "Canvas Title", 800, 600);

    TLegend legTauStar = getLegend();
    legTauStar.AddEntry(h_up_dsStarTau.GetPtr(), "#pm 1#sigma envelope", "L");
    legTauStar.AddEntry(h_dsStarTau.GetPtr(), "central value", "L");

    h_up_dsStarTau->GetYaxis()->SetTitle("events");
    h_up_dsStarTau->GetXaxis()->SetTitle("Q^{2}(GeV^{2})");

    TText* textTauStar         = new TLatex(0.2, 0.85, "B_{s} #rightarrow D*_{s} #tau #nu");
    textTauStar->SetNDC();   
    textTauStar->SetTextFont(43);
    textTauStar->SetTextSize(25);

    h_up_dsStarTau->Draw("HIST");
    h_dsStarTau->Draw("HIST SAME");
    h_down_dsStarTau->Draw("HIST SAME");
    legTauStar.Draw("SAME");
    textTauStar->Draw("SAME");

    canvTauStar->Update();
    canvTauStar->SaveAs(toSave_dsStarTau.c_str());
    }
 
}   


