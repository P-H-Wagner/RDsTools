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

using namespace std;

// define input files
string getInputFile(string signal){

  string fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/" + signal + "/*";

  return fin_str;

}

TLegend getLegend(){

  TLegend legend(0.5,0.8,0.9,0.9);
  return legend;

}
int main(int nargs, char* args[]){

  gStyle->SetOptStat(0);

  //variable for which we want to compute the average weight
  string var(args[2]); //f.e. q2

  //load average weights and settings 
  YAML::Node average_weights = YAML::LoadFile("average_weights.yaml");
  YAML::Node settings = YAML::LoadFile("average_models.yaml");
  int bins   = settings[var]["bins"].as<int>();
  double min = settings[var]["xmin"].as<double>();
  double max = settings[var]["xmax"].as<double>();

  cout << getInputFile(args[1]) << endl;

  ROOT::RDataFrame* df = nullptr;

  try{
    df = new ROOT::RDataFrame("tree",  getInputFile(args[1]).c_str());  
  }
  catch(const exception& e){ cout << "no file found" << endl; exit(1); }


  auto df_filt  = df->Filter("gen_sig == 0");
  auto h        = df_filt.Histo1D({"h", "", bins, min, max}, var, "central_w");

  h->Scale(1 / average_weights["central_w"].as<float>());
  h->SetLineColor(kBlack); 

  vector<string> directions = {"up","down"};
  for (size_t i = 0; i < 10; i++){
    string key = "e"+to_string(i)+".pdf";
    const char* toSave = key.c_str();

    auto h_up   = df_filt.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_up"); 
    h_up->Scale(1 / average_weights[ "e"+to_string(i)+"_up"].as<float>()); 
    h_up->SetLineColor(kRed); 

    auto h_down = df_filt.Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_down"); 
    h_down->Scale(1 / average_weights[ "e"+to_string(i)+"_down"].as<float>()); 
    h_down->SetLineColor(kBlue); 

    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    h->Draw("HIST");
    h_up->Draw("HIST SAME");
    h_down->Draw("HIST SAME");
    canvas->Update();
    canvas->SaveAs(toSave);

    }

 
 
  //save it
  std::ofstream fout("average_weights.yaml");
  fout << average_weights;

}


