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

  gStyle->SetOptStat(0);

  ROOT::RDataFrame* df_dsMu_tot     = nullptr;
  ROOT::RDataFrame* df_dsTau_tot    = nullptr;
  ROOT::RDataFrame* df_dsStarMu_tot = nullptr;
  ROOT::RDataFrame* df_dsStarTau_tot= nullptr;

  // load gen production file
  try{

    // cohen paper
    //df_dsMu_tot      = new ROOT::RDataFrame("tree",  getInputFile("dsmu_"      +string(args[1]) + "_13_12_2024_14_17_38").c_str());  
    //df_dsTau_tot     = new ROOT::RDataFrame("tree",  getInputFile("dstau_"     +string(args[1]) + "_13_12_2024_14_24_08").c_str());  
    //df_dsStarMu_tot  = new ROOT::RDataFrame("tree",  getInputFile("dsstarmu_"  +string(args[1]) + "_13_12_2024_14_20_50").c_str());  
    //df_dsStarTau_tot = new ROOT::RDataFrame("tree",  getInputFile("dsstartau_" +string(args[1]) + "_13_12_2024_14_24_19").c_str());  

    df_dsMu_tot      = new ROOT::RDataFrame("tree",  getInputFile("dsmu_default_27_01_2025_16_10_39")    .c_str()); // bcl first:dsmu_default_24_01_2025_18_31_27 
    df_dsTau_tot     = new ROOT::RDataFrame("tree",  getInputFile("dstau_default_27_01_2025_16_10_47")   .c_str()); // bcl first:dstau_default_24_01_2025_18_31_37
    df_dsStarMu_tot  = new ROOT::RDataFrame("tree",  getInputFile("dsstarmu_BGLVar_31_01_2025_06_34_38") .c_str()); // old: 13_01_2025_20_36_53 #new unbound: 13_01_2025_14_56_07 #16_01_2025_09_59_15 
    df_dsStarTau_tot = new ROOT::RDataFrame("tree",  getInputFile("dsstartau_BGLVar_31_01_2025_06_34_50").c_str()); // old: 3_01_2025_20_36_42 #new unbound: 13_01_2025_14_56_33 #16_01_2025_09_59_25
 
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
  ofstream fout("average_weights.yaml");
  fout << average_weights;

}


