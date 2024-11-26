#include "Hammer/Tools/HammerRoot.hh"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h>
#include <yaml-cpp/yaml.h>

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

  //load bins, start and stop (produced with histModels.py)
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

  //where to save the average weights
  YAML::Node average_weights;

  //auto h                = df->Histo1D({"h", "", bins, min, max}, var, "central_w");
  float average_central = df->Mean("central_w").GetValue(); 

  cout << "Integral of weighted central_w is: " << average_central;
  average_weights["central_w"] = average_central; 

  vector<string> directions = {"up","down"};
  for (size_t i = 0; i < 10; i++){
    for(auto dir:directions){
      //auto h_dummy = df->Histo1D({"h", "", bins, min, max}, var, "e"+to_string(i)+"_"+dir); 
      float average_var =  df->Mean("e"+to_string(i)+"_"+dir).GetValue();//h_dummy.GetPtr()->Integral();
      cout << "\nIntegral of weighted e" << i << "_"<<dir<<" is: \n" << average_var;
      average_weights["e"+to_string(i)+"_"+dir] = average_var;
    }

  }
 
  //save it
  ofstream fout("average_weights.yaml");
  fout << average_weights;

}


