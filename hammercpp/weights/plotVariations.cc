#include "Hammer/Tools/HammerRoot.hh"
#include "TTree.h"
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
// args[1] = BGLVar, ... (i.e. the target we want to plot with weights)            //
// args[2] = q2_coll, pt_miss, ... (i.e. the variable we want to plot with weights)       //
////////////////////////////////////////////////////////////////////////////////////////////


using namespace std;

// define input files
string getInputFile(string signal){

  string fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/signal_" + signal + "/*";

  return fin_str;

}

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}


void plotWeightEffect(auto h_wout, auto h, const char* key){


  TCanvas *c = new TCanvas("canvas", "Canvas Title", 800, 600);

  h->SetLineStyle(7);
 
  h_wout->SetMaximum(h_wout->GetMaximum()*1.4); 
  h_wout->Draw("HIST");
  h->Draw("HIST SAME");

  TLegend leg = getLegend();
  leg.AddEntry(h_wout.GetPtr(), "Unweighted", "L");
  leg.AddEntry(h.GetPtr(), "Weighted","L");
  leg.Draw("SAME");

  c->Update();
  c->SaveAs(key);


}


int main(int nargs, char* args[]){

  gStyle->SetOptStat(0);

  //variable for which we want to compute the average weight
  string var(args[2]); //f.e. q2_coll

  //load average weights and settings 

  YAML::Node average_weights = YAML::LoadFile("average_weights.yaml"); //from getAverageWeights.cc

  YAML::Node settings        = YAML::LoadFile("average_models.yaml"); //from histModels.py
  int bins   = settings[var]["bins"].as<int>();
  double min = settings[var]["xmin"].as<double>();
  double max = settings[var]["xmax"].as<double>();

  ROOT::RDataFrame* df = nullptr;

  try{
    df = new ROOT::RDataFrame("tree",  getInputFile(args[1]).c_str());  
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


  h_dsMu->Scale(     1 / average_weights["central_w_dsmu"].as<float>());
  h_dsTau->Scale(    1 / average_weights["central_w_dstau"].as<float>());
  h_dsStarMu->Scale( 1 / average_weights["central_w_dsstarmu"].as<float>());
  h_dsStarTau->Scale(1 / average_weights["central_w_dsstartau"].as<float>());


  //First, see the effect of the central weight
 
  plotWeightEffect(h_dsMu_wout,      h_dsMu ,      "dsmu_weight_effect.pdf");
  plotWeightEffect(h_dsTau_wout,     h_dsTau ,     "dstau_weight_effect.pdf");
  plotWeightEffect(h_dsStarMu_wout,  h_dsStarMu ,  "dsstarmu_weight_effect.pdf");
  plotWeightEffect(h_dsStarTau_wout, h_dsStarTau , "dsstartau_weight_effect.pdf");


  h_dsMu->SetLineColor(kBlack); 
  h_dsTau->SetLineColor(kBlack); 
  h_dsStarMu->SetLineColor(kBlack); 
  h_dsStarTau->SetLineColor(kBlack); 

  vector<string> directions = {"up","down"};

  for (size_t i = 0; i < 10; i++){

    string key_dsMu      = "e"+to_string(i)+"_dsmu.pdf";
    string key_dsTau     = "e"+to_string(i)+"_dstau.pdf";
    string key_dsStarMu  = "e"+to_string(i)+"_dsstarmu.pdf";
    string key_dsStarTau = "e"+to_string(i)+"_dsstartau.pdf";

    const char* toSave_dsMu      = key_dsMu.c_str();
    const char* toSave_dsTau     = key_dsTau.c_str();
    const char* toSave_dsStarMu  = key_dsStarMu.c_str();
    const char* toSave_dsStarTau = key_dsStarTau.c_str();

    auto h_up_dsMu       = df_dsMu.Histo1D(     {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsTau      = df_dsTau.Histo1D(    {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsStarMu   = df_dsStarMu.Histo1D( {"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    auto h_up_dsStarTau  = df_dsStarTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 

    h_up_dsMu->Scale(      1.0 / average_weights[ "e"+to_string(i)+"_up_dsmu"].as<float>()); 
    h_up_dsTau->Scale(     1.0 / average_weights[ "e"+to_string(i)+"_up_dstau"].as<float>()); 
    h_up_dsStarMu->Scale(  1.0 / average_weights[ "e"+to_string(i)+"_up_dsstarmu"].as<float>()); 
    h_up_dsStarTau->Scale( 1.0 / average_weights[ "e"+to_string(i)+"_up_dsstartau"].as<float>()); 

    h_up_dsMu->SetLineColor(kRed); 
    h_up_dsTau->SetLineColor(kRed); 
    h_up_dsStarMu->SetLineColor(kRed); 
    h_up_dsStarTau->SetLineColor(kRed); 

    auto h_down_dsMu       = df_dsMu.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsTau      = df_dsTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsStarMu   = df_dsStarMu.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    auto h_down_dsStarTau  = df_dsStarTau.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 

    h_down_dsMu->Scale(     1.0 / average_weights[ "e"+to_string(i)+"_down_dsmu"].as<float>()); 
    h_down_dsTau->Scale(    1.0 / average_weights[ "e"+to_string(i)+"_down_dstau"].as<float>()); 
    h_down_dsStarMu->Scale( 1.0 / average_weights[ "e"+to_string(i)+"_down_dsstarmu"].as<float>()); 
    h_down_dsStarTau->Scale(1.0 / average_weights[ "e"+to_string(i)+"_down_dsstartau"].as<float>()); 

    h_down_dsMu->SetLineColor(kBlue); 
    h_down_dsTau->SetLineColor(kBlue); 
    h_down_dsStarMu->SetLineColor(kBlue); 
    h_down_dsStarTau->SetLineColor(kBlue); 

    TCanvas *canvMu = new TCanvas("canvas", "Canvas Title", 800, 600);
    h_up_dsMu->Draw("HIST");
    h_dsMu->Draw("HIST SAME");
    h_down_dsMu->Draw("HIST SAME");
    canvMu->Update();
    canvMu->SaveAs(toSave_dsMu);

    cout << " integral  " <<  h_up_dsMu->Integral();
    cout << " integral  " <<  h_dsMu->Integral();

    /*
    TCanvas *canvTau = new TCanvas("canvas", "Canvas Title", 800, 600);
    h_up_dsTau->Draw("HIST");
    h_dsTau->Draw("HIST SAME");
    h_down_dsTau->Draw("HIST SAME");
    canvTau->Update();
    canvTau->SaveAs(toSave_dsTau);

    TCanvas *canvMuStar = new TCanvas("canvas", "Canvas Title", 800, 600);
    h_up_dsStarMu->Draw("HIST");
    h_dsStarMu->Draw("HIST sAME");
    h_down_dsStarMu->Draw("HIST SAME");
    canvMuStar->Update();
    canvMuStar->SaveAs(toSave_dsStarMu);


    TCanvas *canvTauStar = new TCanvas("canvas", "Canvas Title", 800, 600);
    h_up_dsStarTau->Draw("HIST");
    h_dsStarTau->Draw("HIST SAME");
    h_down_dsStarTau->Draw("HIST SAME");
    canvTauStar->Update();
    canvTauStar->SaveAs(toSave_dsStarTau);
    */
    }
 
}   


