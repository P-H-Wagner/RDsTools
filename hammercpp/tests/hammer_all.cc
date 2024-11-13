//HAMMER and HAMMER ROOT libs 
#include "Hammer/Tools/HammerRoot.hh"
#include "TFile.h"
#include "TTree.h"
#include "Hammer/Hammer.hh"
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "Hammer/Particle.hh"
#include "Hammer/IndexTypes.hh"
#include "Hammer/Math/FourMomentum.hh"
#include "Hammer/Process.hh"

#include <yaml-cpp/yaml.h>
#include <boost/math/special_functions/sign.hpp>
//#include "Math/Point3D.h"

#include <iostream>
using namespace std;


class particle{

  public:
    double pt;
    double eta;
    double phi;
    double mass;
    int    pdgid;

};

void printVariations(const std::map<std::string, double>& variations) {
    for (const auto& entry : variations) {
        std::cout << entry.first << " : " << entry.second << std::endl;
    }
}

// get BGL parameter names
vector<string> getBGLParameters(string channel){

  vector<string> pars;

  if (channel == "BsDs") pars = {"a00","a01","a02","a03","ap0","ap1","ap2","ap3"}; //scalar
  else                   pars = {"a0","a1","a2","b0","b1","b2","c1","c2", "d0", "d1"}; //vectorial

  return pars; 
}

// get BGL variations
string getBGLVariations(string channel){

  string path;
  if (channel == "BsDs") path = "/work/pahwagne/RDsTools/hammer/bgl_scalar_variations.yaml";
  else                   path = "/work/pahwagne/RDsTools/hammer/bgl_vector_variations.yaml";

  return path;
}

// input files

const char* getInputFile (const char* signal){

  // old
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/29_08_2024_09_46_52/all_signals_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/09_10_2024_14_41_07/dsstartau_flatChunk_0.root"; 
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/08_10_2024_08_30_13/dstau_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_10_2024_14_06_09/dsmu_flatChunk_0.root";

  const char* fin;

  // use string comparison to compare char* 
  if      (strcmp(signal, "dsmu") == 0)      fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_18_37_40/dsmu_flatChunk_1.root";
  else if (strcmp(signal, "dsstarmu") == 0)  fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_21_02_57/dsstarmu_flatChunk_1.root";
  else if (strcmp(signal, "dstau") == 0)     fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_21_06_31/dstau_flatChunk_1.root";
  else if (strcmp(signal, "dsstartau") == 0) fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_22_07_42/dsstartau_flatChunk_1.root";
  else{
    cout << "Not a valid signal channel, choose dsmu, dsstarmu, dstau or dsstartau. Aborting ...." << endl;  
    exit(1);
  }

  return fin;
}


// prepare Hammer object
void prepareHammer(Hammer::Hammer& ham, string decay, string hadronicDecay, string inputScheme, string targetScheme){

  //set decay
  ham.includeDecay(decay);
  
  //set input scheme
  cout << hadronicDecay << " " << inputScheme << " " << endl;
  ham.setFFInputScheme({{hadronicDecay, inputScheme}});

  //define target
  cout << "scheme= " + decay << " " << hadronicDecay << " " << targetScheme << endl;
  ham.addFFScheme( decay , {{hadronicDecay, targetScheme}});

  //define units and start run
  ham.setUnits("GeV");
  ham.initRun();


}

void defineDecay(Hammer::Hammer& ham, Hammer::Process& pr, particle b, particle l, particle h){

  ROOT::Math::PtEtaPhiMVector bTlv (b.pt, b.eta, b.phi,  b.mass);
  ROOT::Math::PtEtaPhiMVector lTlv (l.pt, l.eta, l.phi, l.mass);
  ROOT::Math::PtEtaPhiMVector hTlv (h.pt, h.eta, h.phi, h.mass);
  ROOT::Math::PtEtaPhiMVector nuTlv = bTlv - lTlv - hTlv;

  //set pdg id of neutrino:
  int nu_pdgid = (abs(l.pdgid) + 1) * (-1) * boost::math::sign(l.pdgid);
  cout << int(b.pdgid) << " " << int(l.pdgid) << " " << int(h.pdgid) << " " << int(nu_pdgid) << endl;

  Hammer::Particle bPar(  Hammer::FourMomentum( bTlv.e(),  bTlv.px(),  bTlv.py(),  bTlv.pz()),  b.pdgid) ;
  Hammer::Particle lPar(  Hammer::FourMomentum( lTlv.e(),  lTlv.px(),  lTlv.py(),  lTlv.pz()),  l.pdgid) ;
  Hammer::Particle hPar(  Hammer::FourMomentum( hTlv.e(),  hTlv.px(),  hTlv.py(),  hTlv.pz()),  h.pdgid) ;
  Hammer::Particle nuPar( Hammer::FourMomentum( nuTlv.e(), nuTlv.px(), nuTlv.py(), nuTlv.pz()), nu_pdgid);

  // get all involved indices    
  Hammer::ParticleIndex bIdx  = pr.addParticle(bPar);
  Hammer::ParticleIndex lIdx  = pr.addParticle(lPar);
  Hammer::ParticleIndex hIdx  = pr.addParticle(hPar);
  Hammer::ParticleIndex nuIdx = pr.addParticle(nuPar);
  
  // daughters
  Hammer::ParticleIndices daus = {hIdx, lIdx, nuIdx};

  // add decay vertex
  pr.addVertex(bIdx, daus);

  // add process
  auto processId = ham.addProcess(pr);

  // process event
  ham.processEvent();

};

double getCentralWeight(Hammer::Hammer& ham, vector<string> pars, string decay, string hadronicDecay, string targetScheme){

  map<string, double> variations_central;

  // Set all variations to zero to first calculate the central weights
  for(auto name: pars) variations_central["delta_" + name] = 0;

  // Set these variations
  ham.setFFEigenvectors(hadronicDecay, targetScheme, variations_central);

  // print central weight
  double weight = ham.getWeight(decay);

  cout << "central weight " << weight << endl;

  return weight;


};

vector<double> getVariationWeights(Hammer::Hammer& ham, vector<string> pars, YAML::Node vars, string decay, string hadronicDecay, string targetScheme){

  map<string, double> variations;
  vector<string> directions = {"up","down"};
  vector<double> weights;

  //Now we loop over every parameter, and calculate its up/down variation weight
  for(auto name: pars){
    for(auto dir: directions){

      cout << "setting weights for parameter " << name << "in direction " << dir << endl;

      //now for parameter <name> in direction <dir> we set the variations:
      for(auto varPar: pars){

        variations["delta_"+varPar] = vars[name][dir]["delta_" +varPar].as<double>();
        cout << "delta_"+varPar+" has value " << vars[name][dir]["delta_" +varPar].as<double>() <<endl;

      }

      printVariations(variations);

      ham.setFFEigenvectors(hadronicDecay, targetScheme, variations);
      double weight = ham.getWeight(decay);
      cout << " ====> Gives variation weight delta_" +name+dir + " " << weight <<endl;
  
      weights.push_back(weight);

    }


  }

  return weights;

};

int main(int nargs, char* args[]){

  // check if command line arg is given
  if (nargs < 2) {
    //./hammer_all.cc counts already as 1!
    cout << "Please specify signal: dsmu, dstau, dsstarmu, dsstartau";
    exit(1);
  }



  // start setting up hammer
  Hammer::IOBuffer fbBuffer;
  
  // need one hammer instance for every signal!
  Hammer::Hammer hamDsTau;
  Hammer::Hammer hamDsStarTau;
  Hammer::Hammer hamDsMu;
  Hammer::Hammer hamDsStarMu;
  
  // define input/output FF for every signal
  prepareHammer(hamDsTau,     "BsDsTauNu",  "BsDs",  "ISGW2", "BGLVar");
  prepareHammer(hamDsStarTau, "BsDs*TauNu", "BsDs*", "ISGW2", "BGLVar");
  prepareHammer(hamDsMu,      "BsDsMuNu",   "BsDs",  "CLN",   "BGLVar");
  prepareHammer(hamDsStarMu,  "BsDs*MuNu",  "BsDs*", "CLN",   "BGLVar");
  
  //get parameter and variations
  vector<string> parsDsTau     = getBGLParameters("BsDs");
  vector<string> parsDsStarTau = getBGLParameters("BsDs*");
  vector<string> parsDsMu      = getBGLParameters("BsDs");
  vector<string> parsDsStarMu  = getBGLParameters("BsDs*");
  
  string pathDsTau             = getBGLVariations("BsDs"); 
  string pathDsStarTau         = getBGLVariations("BsDs*"); 
  string pathDsMu              = getBGLVariations("BsDs"); 
  string pathDsStarMu          = getBGLVariations("BsDs*"); 
  
  YAML::Node varsDsTau         = YAML::LoadFile(pathDsTau);
  YAML::Node varsDsStarTau     = YAML::LoadFile(pathDsStarTau);
  YAML::Node varsDsMu          = YAML::LoadFile(pathDsMu);
  YAML::Node varsDsStarMu      = YAML::LoadFile(pathDsStarMu);
  
  //load root File for signal given in command line (dsmu, dstau, dsstar mu, dsstar tau)
  unique_ptr<TFile> fin(TFile::Open( getInputFile(args[1]) ));
  
  if (!fin || fin->IsZombie()) {
    cout << "Could not open file!" << endl;
    exit(-1);
  }
  
  //get root tree
  unique_ptr<TTree> tree(fin->Get<TTree>("tree"));
  if (!tree){
    cout << "Could not open tree!" << endl;
    exit(-1);
  }
  
  
  //max nr of events to loop over
  int maxevents = 1;
  
  //get pt, eta, phi of involved particles
  double bs_pt          = 0;
  double bs_eta         = 0;
  double bs_phi         = 0;
  double bs_m           = 0;
  double bs_pdgid       = 0;

  double mu_pt          = 0;
  double mu_eta         = 0;
  double mu_phi         = 0;
  double mu_m           = 0;
  double mu_pdgid       = 0; 
  
  double tau_pt         = 0;
  double tau_eta        = 0;
  double tau_phi        = 0;
  double tau_m          = 0;
  double tau_pdgid      = 0; 
  
  double ds_pt          = 0;
  double ds_eta         = 0;
  double ds_phi         = 0;
  double ds_m           = 0;
  double ds_charge      = 0;
  double ds_pdgid       = 0;
  
  double dsStar_pt      = 0;
  double dsStar_eta     = 0;
  double dsStar_phi     = 0;
  double dsStar_m       = 0;
  double dsStar_pdgid   = 0;
  
  double nu_pdgid       = 0; 
  
  double sig            = 0;
  double event          = 0;
  
  tree->SetBranchAddress("gen_bs_pt",        &bs_pt);
  tree->SetBranchAddress("gen_bs_eta",       &bs_eta);
  tree->SetBranchAddress("gen_bs_phi",       &bs_phi);
  tree->SetBranchAddress("gen_bs_m",         &bs_m);
  tree->SetBranchAddress("gen_bs_pdgid",     &bs_pdgid);

  tree->SetBranchAddress("gen_mu_pt",        &mu_pt);
  tree->SetBranchAddress("gen_mu_eta",       &mu_eta);
  tree->SetBranchAddress("gen_mu_phi",       &mu_phi);
  tree->SetBranchAddress("gen_mu_m",         &mu_m);
  tree->SetBranchAddress("gen_mu_pdgid",     &mu_pdgid);
  
  tree->SetBranchAddress("gen_tau_pt",       &tau_pt);
  tree->SetBranchAddress("gen_tau_eta",      &tau_eta);
  tree->SetBranchAddress("gen_tau_phi",      &tau_phi);
  tree->SetBranchAddress("gen_tau_m",        &tau_m);
  tree->SetBranchAddress("gen_tau_pdgid",    &tau_pdgid);
  
  tree->SetBranchAddress("gen_ds_pt",        &ds_pt);
  tree->SetBranchAddress("gen_ds_eta",       &ds_eta);
  tree->SetBranchAddress("gen_ds_phi",       &ds_phi);
  tree->SetBranchAddress("gen_ds_m",         &ds_m);
  tree->SetBranchAddress("gen_ds_pdgid",     &ds_pdgid);
  tree->SetBranchAddress("gen_ds_charge",    &ds_charge);
  
  tree->SetBranchAddress("gen_dsStar_pt",    &dsStar_pt);
  tree->SetBranchAddress("gen_dsStar_eta",   &dsStar_eta);
  tree->SetBranchAddress("gen_dsStar_phi",   &dsStar_phi);
  tree->SetBranchAddress("gen_dsStar_m",     &dsStar_m);
  tree->SetBranchAddress("gen_dsStar_pdgid", &dsStar_pdgid);
  
  tree->SetBranchAddress("gen_sig",          &sig);
  tree->SetBranchAddress("event",            &event);
  
  cout << "To be processed events: " << maxevents << endl;
  
  for(size_t i = 0; i < maxevents; ++i){
  
    if(i%1000 == 0) cout << " ====> Processing event " << i << endl;
  
    //start event
    hamDsTau.initEvent();
    hamDsStarTau.initEvent();
  
    // load entry i into memory 
    tree->GetEntry(i);
  
    cout << "NEW EVENT with id " << int(sig) << endl;
    cout << " bs pt " << bs_pt << " tau pt " << tau_pt << endl;
  
    // Define the decays for which you set up hammer    
    set<int> validDecays = {0,1,10,11};
  
    if (validDecays.find(int(sig)) == validDecays.end()){
      cout << "Signal ID is " << sig << ", skipping event ..." << endl;
      continue;
    }
  
    //Define a class for each particle
    //All signals are of type Bs -> hc + lep + neutrino
    particle bs;
    particle lep;
    particle hc;
  
    //Bs is the same for all signals
    ROOT::Math::PtEtaPhiMVector bsTlv (bs_pt,  bs_eta,  bs_phi,  bs_m);
    Hammer::Particle bsPar( Hammer::FourMomentum(bsTlv.e(),  bsTlv.px(),  bsTlv.py(),  bsTlv.pz()),  int(bs_pdgid));
    bs.pt  = bs_pt;     bs.phi  = bs_phi;     bs.eta  = bs_eta;     bs.mass  = bs_m;     bs.pdgid  = bs_pdgid;
  
    //create the process
    Hammer::Process process;
 
    if (sig == 0){
      //////////////////////////
      // Bs -> Ds + Mu + Nu   //
      //////////////////////////

      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsMu, process, bs, lep, hc);
  
      // get central weight
      double central_weight = getCentralWeight(hamDsMu, parsDsMu, "BsDsMuNu", "BstoDs", "BGLVar");
     
      //get variation weights
      vector<double> weights = getVariationWeights(hamDsMu, parsDsMu, varsDsMu, "BsDsMuNu", "BstoDs", "BGLVar");
  
    }
 
    if (sig == 1){
      //////////////////////////
      // Bs -> Ds + Tau + Nu  //
      //////////////////////////
  
      lep.pt = tau_pt;    lep.phi = tau_phi;    lep.eta = tau_eta;    lep.mass = tau_m;    lep.pdgid = tau_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsTau, process, bs, lep, hc);
  
      // get central weight
      double central_weight = getCentralWeight(hamDsTau, parsDsTau, "BsDsTauNu", "BstoDs", "BGLVar");
     
      //get variation weights
      vector<double> weights = getVariationWeights(hamDsTau, parsDsTau, varsDsTau, "BsDsTauNu", "BstoDs", "BGLVar");
  
    }
 
    if (sig == 10){
      //////////////////////////
      // Bs -> Ds* + Mu + Nu  //
      //////////////////////////
  
      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = dsStar_pt; hc.phi  = dsStar_phi; hc.eta  = dsStar_eta; hc.mass  = dsStar_m; hc.pdgid  = dsStar_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarMu, process, bs, lep, hc);
  
      // get central weight
      double central_weight = getCentralWeight(hamDsStarMu, parsDsStarMu, "BsDs*MuNu", "BstoDs*", "BGLVar");
     
      //get variation weights
      vector<double> weights = getVariationWeights(hamDsStarMu, parsDsStarMu, varsDsStarMu, "BsDs*MuNu", "BstoDs*", "BGLVar");
  
    }
 
    if (sig == 11){
      //////////////////////////
      // Bs -> Ds* + Tau + Nu //
      //////////////////////////
  
      lep.pt = tau_pt;    lep.phi = tau_phi;    lep.eta = tau_eta;    lep.mass = tau_m;    lep.pdgid = tau_pdgid;
      hc.pt  = dsStar_pt; hc.phi  = dsStar_phi; hc.eta  = dsStar_eta; hc.mass  = dsStar_m; hc.pdgid  = dsStar_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarTau, process, bs, lep, hc);
  
      // get central weight
      double central_weight = getCentralWeight(hamDsStarTau, parsDsStarTau, "BsDs*TauNu", "BstoDs*", "BGLVar");
     
      //get variation weights
      vector<double> weights = getVariationWeights(hamDsStarTau, parsDsStarTau, varsDsStarTau, "BsDs*TauNu", "BstoDs*", "BGLVar");
  
    } //closing if
  
  } // closing event loop

} //closing main

