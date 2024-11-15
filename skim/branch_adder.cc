#include <iostream>
#include <ROOT/RDataFrame.hxx>

using namespace std;

//adding branches w.r.t to score trees


/////////////////////////////////////////
// Input single file
const char* getInputFile(TString channel){

  const char* fin;
  if (channel == "sig") fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_26Sep2024_07h46m21s_cons/sig_26Sep2024_07h46m21s_flatChunk_0.root";
  else { cout << "No channel named: "<< channel << endl; exit(1);}

  return fin;

}

/////////////////////////////////////////
// Input chain 
TChain* getInputChain(TString channel){

  const char* fin;
  if (channel == "sig") fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_26Sep2024_07h46m21s_cons/*";
  else { cout << "No channel named: "<< channel << endl; exit(1);}

  TChain* chain = new TChain("tree");
  chain->Add(fin);
  return chain;

}

/////////////////////////////////////////
// Define perp. Ds momentum w.r.t Bs flight direction 
// ugly long function but i dont see a better option whenusing it together with RDataFrame's define method ...
double ds_perp(double ds_pt, double ds_eta, double ds_phi, double ds_m, double pv_x, double pv_y, double pv_z, double sv_x, double sv_y, double sv_z ){

  TLorentzVector dsTlv;
  dsTlv.SetPtEtaPhiM(ds_pt, ds_eta, ds_phi, ds_m);
  TVector3 bFlightDir;
  bFlightDir.SetXYZ(pv_x - sv_x, pv_y - sv_y, pv_z - sv_z);
 
  double dsPerp   = dsTlv.Vect().Perp(bFlightDir);
  
  return dsPerp; 

}

/////////////////////////////////////////
// Define corrected mass 
// ugly long function but i dont see a better option whenusing it together with RDataFrame's define method ...
double m_corr(double dsMu_m, double ds_pt, double ds_eta, double ds_phi, double ds_m, double pv_x, double pv_y, double pv_z, double sv_x, double sv_y, double sv_z ){

  TLorentzVector dsTlv;
  dsTlv.SetPtEtaPhiM(ds_pt, ds_eta, ds_phi, ds_m);
  TVector3 bFlightDir;
  bFlightDir.SetXYZ(pv_x - sv_x, pv_y - sv_y, pv_z - sv_z);
 
  double dsPerp   = dsTlv.Vect().Perp(bFlightDir);
  double mCorr = std::sqrt(std::pow(dsMu_m,2) + std::pow(dsPerp,2)) + dsPerp; 

  return mCorr; 

}


/////////////////////////////////////////
// Main func 
void branch_adder(TString channel){
  cout << "hi" << endl;


  //get input chain
  TChain* chain = getInputChain(channel);

  int maxevents = chain->GetEntries();

  //load into df
  ROOT::RDataFrame df_full(*chain);

  // get input file 
  //const char* fin = getInputFile(channel);
  //load it into df
  //ROOT::RDataFrame df_full("tree", fin);


  //load only part of df to debug
  auto df = df_full.Range(0,maxevents);
  auto df_with_columns = df.Define("test", "mu_pt*2")
  .Define("phiPi_perp",ds_perp, {"phiPi_pt","phiPi_eta","phiPi_phi","phiPi_m","pv_x","pv_y","pv_z","sv_x","sv_y","sv_z" })
  .Define("m_corr", m_corr, {"dsMu_m","phiPi_pt","phiPi_eta","phiPi_phi","phiPi_m","pv_x","pv_y","pv_z","sv_x","sv_y","sv_z" });

  df_with_columns.Snapshot("tree", "/scratch/pahwagne/adder/test.root");
}
