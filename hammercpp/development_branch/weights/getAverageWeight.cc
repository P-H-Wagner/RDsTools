#include "Hammer/Tools/HammerRoot.hh"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include <ROOT/RDataFrame.hxx>
#include <TStyle.h>
#include <yaml-cpp/yaml.h>

/////////////////////////////////////////////////////////////////////////////
// Command line args:                                                      //
/////////////////////////////////////////////////////////////////////////////

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

  cout << args[0] << endl; 
  cout << args[1] << endl; //dsmu
  cout << args[2] << endl; //dstau
  cout << args[3] << endl; //dsstarmu
  cout << args[4] << endl; //dsstartau
  cout << args[5] << endl; //datetime

  gStyle->SetOptStat(0);

  //load hammer trees
  string dsmu_fin       = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/" + string(args[1]) + "/*";
  string dstau_fin      = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/" + string(args[2]) + "/*";
  string dsstarmu_fin   = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/" + string(args[3]) + "/*";
  string dsstartau_fin  = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/hammer/25/" + string(args[4]) + "/*";


  ROOT::RDataFrame* df_dsMu_tot     = nullptr;
  ROOT::RDataFrame* df_dsTau_tot    = nullptr;
  ROOT::RDataFrame* df_dsStarMu_tot = nullptr;
  ROOT::RDataFrame* df_dsStarTau_tot= nullptr;

  // load gen production file
  try{

    df_dsMu_tot      = new ROOT::RDataFrame("tree", dsmu_fin      ); 
    df_dsTau_tot     = new ROOT::RDataFrame("tree", dstau_fin     ); 
    df_dsStarMu_tot  = new ROOT::RDataFrame("tree", dsstarmu_fin  ); 
    df_dsStarTau_tot = new ROOT::RDataFrame("tree", dsstartau_fin ); 
 
  }
  catch(const exception& e){ cout << "no file found" << endl; exit(1); }

  auto df_dsMu       = df_dsMu_tot     ->Filter("gen_sig == 0").Filter([](float x) { return !std::isnan(x); }, {"central_w"});
  auto df_dsTau      = df_dsTau_tot    ->Filter("gen_sig == 1").Filter([](float x) { return !std::isnan(x); }, {"central_w"});
  auto df_dsStarMu   = df_dsStarMu_tot ->Filter("gen_sig == 10").Filter([](float  x) { return !std::isnan(x); }, {"central_w"});
  auto df_dsStarTau  = df_dsStarTau_tot->Filter("gen_sig == 11").Filter([](float  x) { return !std::isnan(x); }, {"central_w"});

  //where to save the average weights
  YAML::Node average_weights;

  // access central weights and take average
  double average_central_dsmu      = df_dsMu.Mean("central_w").GetValue(); 
  double average_central_dstau     = df_dsTau.Mean("central_w").GetValue(); 
  double average_central_dsstarmu  = df_dsStarMu.Mean("central_w").GetValue(); 
  double average_central_dsstartau = df_dsStarTau.Mean("central_w").GetValue(); 

  cout << "Average of dsmu weighted central_w is: "      << average_central_dsmu;
  cout << "Average of dstau weighted central_w is: "     << average_central_dstau;
  cout << "Average of dsstarmu weighted central_w is: "  << average_central_dsstarmu;
  cout << "Average of dsstartau weighted central_w is: " << average_central_dsstartau;

  average_weights["central_w_dsmu"]      = average_central_dsmu; 
  average_weights["central_w_dstau"]     = average_central_dstau; 
  average_weights["central_w_dsstarmu"]  = average_central_dsstarmu; 
  average_weights["central_w_dsstartau"] = average_central_dsstartau; 

  // access variation weights and take average
  vector<string> directions = {"up","down"};

  for (size_t i = 1; i < 11; i++){
    for(auto dir:directions){
      float average_var_dsMu       =  df_dsMu.Mean("e"+to_string(i)+"_"+dir).GetValue();
      float average_var_dsTau      =  df_dsTau.Mean("e"+to_string(i)+"_"+dir).GetValue();
      float average_var_dsStarMu   =  df_dsStarMu.Mean("e"+to_string(i)+"_"+dir).GetValue();
      float average_var_dsStarTau  =  df_dsStarTau.Mean("e"+to_string(i)+"_"+dir).GetValue();

      cout << "\nAverage of dsmu variation weight " << i << "_"<<dir<<" is: \n" << average_var_dsMu;
      cout << "\nAverage of dstau variation weight " << i << "_"<<dir<<" is: \n" << average_var_dsTau;
      cout << "\nAverage of dsstarmu variation weight " << i << "_"<<dir<<" is: \n" << average_var_dsStarMu;
      cout << "\nAverage of dsstartau variation weight " << i << "_"<<dir<<" is: \n" << average_var_dsStarTau;


      average_weights["e"+to_string(i)+"_"+dir+"_dsmu"] = average_var_dsMu;
      average_weights["e"+to_string(i)+"_"+dir+"_dstau"] = average_var_dsTau;
      average_weights["e"+to_string(i)+"_"+dir+"_dsstarmu"] = average_var_dsStarMu;
      average_weights["e"+to_string(i)+"_"+dir+"_dsstartau"] = average_var_dsStarTau;

    }
  }
 
  //save it
  string saveYaml = string(args[5]) + "/average_weights.yaml";
  ofstream fout(saveYaml.c_str());
  fout << average_weights;

}


