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

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

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

string getCLNParameters(string decay){

  string paras = "{";

  if (decay == "BsDsMuNu"){

    cout << "Adapting scalar HQET --> CLN" << endl;

    // EvtGen settings
    double evtrho2  = 1.17;
    double evtv1_1  = 1.074;

    // TODO
    double d_m      = 1.96834;
    double b_m      = 5.36688;

    double rC       = d_m / b_m ; 
    double clnRhoSq = evtrho2;
    double clnG1    = evtv1_1 * 2 * sqrt(rC) / (1 + rC);
   
    paras += "RhoSq: " + to_string(clnRhoSq) + ", ";
    paras += "G1: "    + to_string(clnG1);
    paras += "}";    

  }

  else if (decay == "BsDs*MuNu"){

    cout << "Adapting vectorial HQET --> CLN" << endl;

    // EvtGen settings
    double evtrho2   = 1.16;
    double evtha1_1  = 0.921;
    double evtr1     = 1.37;
    double evtr2     = 0.845;
 
    double clnRhoSq  = evtrho2;
    double clnF1     = evtha1_1;
    double clnR1     = evtr1;
    double clnR2     = evtr2;

    paras += "RhoSq: " + to_string(clnRhoSq) + ", ";
    paras += "F1: "    + to_string(clnF1)    + ", ";
    paras += "R1: "    + to_string(clnR1)    + ", ";
    paras += "R2: "    + to_string(clnR2);
    paras += "}";    


  }

  else{

    cout << "Not a valid decay for getCLNParameters(), choose BsDsMuNu or BsDs*MuNu, aborting ..." << endl;
    exit(1);

  }

  return paras;

};

// input files

const char* getInputFile(const char* signal){

  // old
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/29_08_2024_09_46_52/all_signals_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/09_10_2024_14_41_07/dsstartau_flatChunk_0.root"; 
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/08_10_2024_08_30_13/dstau_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_10_2024_14_06_09/dsmu_flatChunk_0.root";

  const char* fin;

  // use string comparison to compare char* 
  if      (strcmp(signal, "dsmu") == 0) {
    //fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_18_37_40/*";
    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/*"; //more stats!

  }
  else if (strcmp(signal, "dsmu_isgw2") == 0) {
    //fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/13_11_2024_07_15_45/*";
    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/17_11_2024_17_26_37/*"; //more stats!
  }

  else if (strcmp(signal, "signal") == 0){

    fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/skimmed/26_07_2024_14_46_03/*";

  }
  //else if (strcmp(signal, "dsstarmu") == 0)  fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_21_02_57/*";
  //else if (strcmp(signal, "dstau") == 0)     fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_21_06_31/*";
  //else if (strcmp(signal, "dsstartau") == 0) fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_22_07_42/*";
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
  ham.addFFScheme("schemeInput", {{hadronicDecay, inputScheme}});
  ham.setFFInputScheme({{hadronicDecay, inputScheme}});

  if ((inputScheme == "CLN") || (targetScheme == "CLN")){

    // if input scheme is CLN, we have to adapt parameters from HQET (EvtGen) -> CLN (Hammer)
    string paras = getCLNParameters(decay);

    // string is of type ham.setOptions("BstoDsBGL: {ChiT: 0.01, ChiL: 0.002}");
    string key;

    if(decay.find("Ds*") != string::npos) key = "BstoDs*CLN: " + paras;
    else                                  key = "BstoDsCLN: "  + paras;
    cout << "the key is: " << key << endl;
    ham.setOptions(key);

  }

  //define target
  cout << "scheme= " + decay << " " << hadronicDecay << " " << targetScheme << endl;
  ham.addFFScheme( "scheme"+decay , {{hadronicDecay, targetScheme}});

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
  //cout << int(b.pdgid) << " " << int(l.pdgid) << " " << int(h.pdgid) << " " << int(nu_pdgid) << endl;
  //cout << "lepton pt: " << l.pt << endl;

  Hammer::Particle bPar(  { bTlv.e(),  bTlv.px(),  bTlv.py(),  bTlv.pz()},  b.pdgid) ;
  Hammer::Particle lPar(  { lTlv.e(),  lTlv.px(),  lTlv.py(),  lTlv.pz()},  l.pdgid) ;
  Hammer::Particle hPar(  { hTlv.e(),  hTlv.px(),  hTlv.py(),  hTlv.pz()},  h.pdgid) ;
  Hammer::Particle nuPar( { nuTlv.e(), nuTlv.px(), nuTlv.py(), nuTlv.pz()}, nu_pdgid);

  // get all involved indices    
  Hammer::ParticleIndex bIdx  = pr.addParticle(bPar);
  Hammer::ParticleIndex lIdx  = pr.addParticle(lPar);
  Hammer::ParticleIndex hIdx  = pr.addParticle(hPar);
  Hammer::ParticleIndex nuIdx = pr.addParticle(nuPar);
  
  // daughters
  Hammer::ParticleIndices daus = {lIdx, hIdx, nuIdx};

  // add decay vertex
  pr.addVertex(bIdx, daus);

  // add process
  auto processId = ham.addProcess(pr);


};

double getCentralWeight(Hammer::Hammer& ham, vector<string> pars, string decay, string hadronicDecay, string targetScheme){

  map<string, double> variations_central;

  // Set all variations to zero to first calculate the central weights
  //for(auto name: pars) variations_central["delta_" + name] = 0.0;

  // Set these variations
  //ham.setFFEigenvectors(hadronicDecay, targetScheme, variations_central);

  // process event
  ham.processEvent();

  // print central weight
  double weight = ham.getWeight("scheme"+decay);

  //cout << "central weight " << weight << endl;

  return weight;


};

map<string, double> getVariationWeights(Hammer::Hammer& ham, vector<string> pars, YAML::Node vars, string decay, string hadronicDecay, string targetScheme){

  map<string, double> variations;
  vector<string> directions = {"up","down"};
  map<string, double> weights;

  int counter = 0;
  //Now we loop over every parameter, and calculate its up/down variation weight
  for(auto name: pars){
    for(auto dir: directions){

      //cout << "setting weights for parameter " << name << "in direction " << dir << endl;

      //now for parameter <name> in direction <dir> we set the variations:
      for(auto varPar: pars){

        variations["delta_"+varPar] = vars[name][dir]["delta_" +varPar].as<double>();
        //cout << "delta_"+varPar+" has value " << vars[name][dir]["delta_" +varPar].as<double>() <<endl;

      }

      //printVariations(variations);

      ham.setFFEigenvectors(hadronicDecay, targetScheme, variations);

      // process event
      ham.processEvent();

      double weight = ham.getWeight(decay);
      //cout << " ====> Gives variation weight delta_" +name+dir + " " << weight <<endl;
 
      //cout << "key is " << "e"+to_string(counter)+"_"+dir << endl; 
      weights["e"+to_string(counter)+"_"+dir] = weight;

    }

    ++counter;

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
  Hammer::Hammer hamDsMu;
  
  // define input/output FF for every signal
  prepareHammer(hamDsMu,      "BsDsMuNu",   "BsDs",  args[2],   args[3]);
  
  //get parameter and variations
  vector<string> parsDsMu      = getBGLParameters("BsDs");
  
  string pathDsMu              = getBGLVariations("BsDs"); 
  
  YAML::Node varsDsMu          = YAML::LoadFile(pathDsMu);
  
  //load root File for signal given in command line (dsmu, dstau, dsstar mu, dsstar tau)
  TChain* tree = new TChain("tree");
  tree->Add( getInputFile(args[1]));
  cout << "tree has entries " << tree->GetEntries() << endl;

  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_0.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_1.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_2.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_3.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_4.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_5.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_6.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_7.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_8.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_9.root");
  //tree->Add("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_10.root");
  /*
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
  */
  //max nr of events to loop over
  int maxevents = tree->GetEntries();
  if (maxevents > tree->GetEntries()) maxevents = tree->GetEntries();

  //define a place holder for the central_weights
  map<int, double> central_weights;
  //and one for the variational weights (nested map, since we have several parameters)
  map<int, map<string,double>> var_weights;

  ROOT::RDataFrame df_full("tree",  getInputFile(args[1]) );
  std::vector<std::string> files = {
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_0.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_1.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_2.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_3.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_4.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_5.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_6.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_7.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_8.root",
    "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/dsmu_flatChunk_9.root"
  }; 

  //ROOT::RDataFrame df_full("tree",  files );
  cout << "full df has entries " << df_full.Count().GetValue() << endl; 

  auto df = df_full.Range(0,maxevents).Filter("gen_sig == 0");
  cout << "filtered df has entries " << df.Count().GetValue() << endl; 
  
  
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
 
  string dest_str = "/scratch/pahwagne/hammer/" + string(args[1]) + "_circ.root";  
  const char* dest = dest_str.c_str();

  TFile outputFile(dest, "RECREATE");
  TTree* outputTree = tree->CloneTree(0);
  float central_w;
  outputTree->Branch("central_w", &central_w, "central_w/F");
 
  for(size_t i = 0; i < maxevents; ++i){
  
    if(i%1000 == 0) cout << " ====> Processing event " << i << endl;
  
    //start event
    hamDsMu.initEvent();
  
    // load entry i into memory 
    tree->GetEntry(i);
  
    //cout << " --------- NEW EVENT of signal " << int(sig) << " -------------" << endl;
    //cout << " bs pt " << bs_pt << " mu pt " << mu_pt << endl;
  
    // Define the decays for which you set up hammer    
    set<int> validDecays = {0};
  
    if (validDecays.find(int(sig)) == validDecays.end()){
      //cout << "Signal ID is " << sig << ", skipping event ..." << event << endl;
      continue;
    }
  
    long long event_int = static_cast<long long>(event);
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

    //placeholders for weights
    double central_weight;
    map<string, double> weights;
 
    if (sig == 0){
      //////////////////////////
      // Bs -> Ds + Mu + Nu   //
      //////////////////////////

      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsMu, process, bs, lep, hc);
  
      // get central weight
      central_weight = getCentralWeight(hamDsMu, parsDsMu, "BsDsMuNu", "BstoDs", args[2]);
     
      //get variation weights
      //weights = getVariationWeights(hamDsMu, parsDsMu, varsDsMu, "BsDsMuNu", "BstoDs", "BGLVar");
  
      //cout << "At index: " << i << " ====> saving " << central_weight << " for event nr " << event_int << endl; 
    }
 
    // assign it to tree variable
    //central_w_leaf = central_weight;
    // fill it
    //cout << "filling now" << endl;
    //treeNew->Fill();

    else {
      central_weight = -1;
      //weight = -1;
    }
    central_w = central_weight;
    outputTree->Fill(); 
    central_weights[event_int] = central_weight; 
    //var_weights[int(event)]     = weights; 

  } // closing event loop

  outputTree->Write();
  outputFile.Close();
 
  //save central weights 
  //df = df.Define("central_w", [central_weights](double nr, double id) { 

  //  long long nr_int = static_cast<long long>(nr);
  //  //cout << "i am at event: " << nr_int << " and  signal ID is: " << int(id) << " and I save " << central_weights.at(nr_int) << endl;
  //  return central_weights.at(nr_int); 
  //}, {"event", "gen_sig"});

  //save variation weights

  // first, get the maximum of the nr of parameters of all channels
  /*
  int n_vars = max({int(parsDsMu.size()), int(parsDsStarMu.size()), int(parsDsTau.size()), int(parsDsStarTau.size())}); 

  
  vector<string> directions = {"up","down"};
  int counter = 0;

  for(size_t i = 0; i < n_vars; ++i){
    for(auto dir:directions){

      string key = "e"+to_string(counter)+"_"+dir;
      //cout << key << endl;

      df = df.Define(key, [counter, key, var_weights](double nr) { 

        if (int(var_weights.at(int(nr)).size() / 2) <= counter ) {
          //cout << "returning -1" << endl; 
          return -1.0;} 
        else { 
          //cout << "accessing weight" << endl; 
          return var_weights.at(int(nr)).at(key);} 
        }, {"event"});

    }
    ++counter;
  }
  */
  //cout << "max number of paras is: " << n_vars << endl;
  cout << "saving to: " << "/scratch/pahwagne/hammer/" + string(args[1]) + "_circ.root" << endl;
  //df.Snapshot("tree", "/scratch/pahwagne/hammer/" + string(args[1]) + "_circ.root");
  //df.Snapshot("tree", string(args[1]) + "_circ.root");


} //closing main

