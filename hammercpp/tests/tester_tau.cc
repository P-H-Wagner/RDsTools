
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
//#include "Math/Point3D.h"

#include <iostream>
using namespace std;

int main(){

  // start setting up hammer
  Hammer::Hammer ham;
  Hammer::IOBuffer fbBuffer;


  ham.includeDecay("BsDsTauNu");
  ham.setFFInputScheme({{"BsDs", "ISGW2"}});
  ham.addFFScheme("Scheme1", {{"BsDs", "BGL"}});
  //ham.addFFScheme("SchemeVar", {{"BsDs", "BGLVar"}});
  //ham.setOptions("BstoDsBGLVar: {a0: [1., 1., 1., 1.], ap: [1., 1., 1., 1.]}");

  ham.setUnits("GeV");
  ham.initRun();

  //load ROOT File
  unique_ptr<TFile> fin(TFile::Open("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/08_10_2024_08_30_13/dstau_flatChunk_0.root"));

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
  int maxevents = 10;

  //get pt, eta, phi of involved particles
  double bs_pt      = 0;
  double bs_eta     = 0;
  double bs_phi     = 0;
  double bs_m       = 5.366;
  double bs_pdgid   = 0;

  double tau_pt     = 0;
  double tau_eta    = 0;
  double tau_phi    = 0;
  double tau_m      = 1.776;
  int    tau_pdgid  = 15; //sign is adapted later

  double ds_pt      = 0;
  double ds_eta     = 0;
  double ds_phi     = 0;
  double ds_m       = 2.112;
  double ds_pdgid   = 0;
  double ds_charge  = 0;

  int    nu_pdgid   = 16; //sign is adapted later

  tree->SetBranchAddress("genSIM_bs_pt",     &bs_pt);
  tree->SetBranchAddress("genSIM_bs_eta",    &bs_eta);
  tree->SetBranchAddress("genSIM_bs_phi",    &bs_phi);
  //tree->SetBranchAddress("genSIM_bs_m",    &bs_m);
  tree->SetBranchAddress("genSIM_bs_pdgid",  &bs_pdgid);
 
  tree->SetBranchAddress("genSIM_tau_pt",    &tau_pt);
  tree->SetBranchAddress("genSIM_tau_eta",   &tau_eta);
  tree->SetBranchAddress("genSIM_tau_phi",   &tau_phi);
  //tree->SetBranchAddress("genSIM_tau_m",   &tau_m);
 
  tree->SetBranchAddress("genSIM_ds_pt",     &ds_pt);
  tree->SetBranchAddress("genSIM_ds_eta",    &ds_eta);
  tree->SetBranchAddress("genSIM_ds_phi",    &ds_phi);
  //tree->SetBranchAddress("genSIM_ds_m",    &ds_m);
  tree->SetBranchAddress("genSIM_ds_pdgid",  &ds_pdgid);
  tree->SetBranchAddress("genSIM_ds_charge", &ds_charge);

  cout << "To be processed events: " << maxevents << endl;

  for(size_t i = 0; i < maxevents; ++i){

    if(i%1000 == 0) cout << "Processing event " << i << endl;

    //start event
    ham.initEvent();

    // load entry i into memory 
    tree->GetEntry(i);

    //create pdg id of tau and nu
    tau_pdgid = 15 * ds_charge;
    nu_pdgid  = 16 * -1 * ds_charge;

    cout << ds_pdgid << " " << bs_pdgid << endl;
    cout << tau_pdgid << " " << nu_pdgid << endl;

    ROOT::Math::PtEtaPhiMVector bsTlv (bs_pt,  bs_eta,  bs_phi,  bs_m); 
    ROOT::Math::PtEtaPhiMVector tauTlv(tau_pt, tau_eta, tau_phi, tau_m); 
    ROOT::Math::PtEtaPhiMVector dsTlv (ds_pt,  ds_eta,  ds_phi,  ds_m); 
    ROOT::Math::PtEtaPhiMVector nuTlv = bsTlv - tauTlv - dsTlv;

    Hammer::Particle bsPar( Hammer::FourMomentum(bsTlv.e(),  bsTlv.px(),  bsTlv.py(),  bsTlv.pz()),  int(bs_pdgid));
    Hammer::Particle tauPar(Hammer::FourMomentum(tauTlv.e(), tauTlv.px(), tauTlv.py(), tauTlv.pz()), int(tau_pdgid));
    Hammer::Particle dsPar( Hammer::FourMomentum(dsTlv.e(),  dsTlv.px(),  dsTlv.py(),  dsTlv.pz()),  int(ds_pdgid));
    Hammer::Particle nuPar( Hammer::FourMomentum(nuTlv.e(),  nuTlv.px(),  nuTlv.py(),  nuTlv.pz()),  int(nu_pdgid));

    //create the process
    Hammer::Process BsDsTauNu;

    // get all involved indices    
    Hammer::ParticleIndex bsIdx  = BsDsTauNu.addParticle(bsPar);
    Hammer::ParticleIndex dsIdx  = BsDsTauNu.addParticle(dsPar);
    Hammer::ParticleIndex tauIdx = BsDsTauNu.addParticle(tauPar);
    Hammer::ParticleIndex nuIdx  = BsDsTauNu.addParticle(nuPar);

    // daughters
    Hammer::ParticleIndices daus = {dsIdx, tauIdx, nuIdx};

    // add decay vertex
    BsDsTauNu.addVertex(bsIdx, daus); 

    // add process
    auto processId = ham.addProcess(BsDsTauNu); 

    // process event
    ham.processEvent();

    // print central weight
    double weights = ham.getWeight("Scheme1");

    // set ff eigenvectors
    // ham.setFFEigenvectors("BsDs", "BGLVar", {{"delta_ap0",1.},{"delta_ap1",1.},{"delta_ap2",1.},{"delta_ap3",1.},{"delta_a00",1.},{"delta_a01",1.},{"delta_a02",1.},{"delta_a03",1.}});

    cout << weights <<endl;


  } 

}

