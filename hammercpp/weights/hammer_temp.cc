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
#include "Hammer/SettingsHandler.hh"

#include <yaml-cpp/yaml.h>
#include <boost/math/special_functions/sign.hpp>
//#include "Math/Point3D.h"

#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

#include <iostream>
#include <string>

#include "bglSettings.h"

#include <sys/stat.h>
#include <unistd.h>

using namespace std;

// before the event loop, lets get the lattice corr and cov matrix!
Eigen::VectorXd latticeCentral = getCentralValues();
auto[latticeCorr, latticeCov]  = getCorrCovMatrix(true);


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

// set BGL parameter names
void setBGLNames(Hammer::Hammer& ham, string channel){

  vector<string> names;
  string key;

  if (channel == "BsDs") {
    names = {"delta_e1","delta_e2","delta_e3","delta_e4","delta_e5","delta_e6","delta_e7", "delta_e8"}; //scalar
    key = "BstoDs"; 
  }
  else{
    names = {"delta_e1","delta_e2","delta_e3","delta_e4","delta_e5","delta_e6","delta_e7", "delta_e8", "delta_e9", "delta_e10"};//vectorial
    key = "BstoDs*";
  }

  ham.renameFFEigenvectors(key, "BGLVar", names);

}

// get BGL parameter names
vector<string> getBGLNames(string channel){

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

// get BGL parameters to set via setOptions()

map<string, map<string, map<string,double>>> getBGLSettings(Hammer::Hammer& ham, string decay,  TLorentzVector d, TLorentzVector b){

  string key = "";
  map<string, map<string, map<string,double>>> variations;
  string paras = "{";

  if ((decay == "BsDsMuNu") or (decay == "BsDsTauNu")){
    //cout << "Setting scalar BGL central values" << endl;
    

    // TODO: import from yaml file
    double a00 = 0.052255946347001495;
    double a01 = -0.16027634967890908;
    double a02 = 0.014141836205563255;
    double a03 = 0.0;
    double ap0 = 0.0017893827864468802;
    double ap1 = -0.004691380424494185;
    double ap2 = -0.015708534616906505;
    double ap3 = 0.0;

    //paras += "a00: " + to_string(a00) + ", ";
    //paras += "a01: " + to_string(a01) + ", ";
    //paras += "a02: " + to_string(a02) + ", ";
    //paras += "a03: " + to_string(a03) + ", ";
    //paras += "ap0: " + to_string(ap0) + ", ";
    //paras += "ap1: " + to_string(ap1) + ", ";
    //paras += "ap2: " + to_string(ap2) + ", ";
    //paras += "ap3: " + to_string(ap3) + ", ";

    paras += "a0: ["  +  to_string(a00)    + ", " + to_string(a01) + ", " + to_string(a02) + ", " + to_string(a03) + "],";
    paras += "ap: ["  +  to_string(ap0)    + ", " + to_string(ap1) + ", " + to_string(ap2) + ", " + to_string(ap3) + "],";

    paras += "}";  

    key = "BstoDsBGLVar: " + paras;
  }

  else{
    //cout << "Setting vector BGL central values" << endl;

    //get q2 and pPrime2
    TLorentzVector q = b - d;
    double q2 = q.M2();
    //cout << "q2 is: " << q2 << endl;
    //cout << "d before: " << d.Pt() << endl;
    boostDsStar(d,b);
    //cout << "d after: " << d.Pt() << endl;
    double pPrime2 = d.Vect().Mag2();   
    //cout << "pPrime2 is: " << pPrime2 <<  endl;
   
    double md      = 1.96834;
    double mb      = 5.36688;
 
    double edp = (mb*mb+md*md-q2)/(2*mb);  
    double momsq = edp*edp-md*md;
    //cout << "momq2 is: " << momsq <<  endl;

    // get change of basis
    Eigen::MatrixXd baseChange     = getBaseChangeMatrix(q2, pPrime2, false);

    //transform the central values!
    Eigen::VectorXd bglCentral  = baseChange * latticeCentral;
    //transform covariance matrix - it is a bilinear form!
    Eigen::MatrixXd bglCov      = baseChange * latticeCov * baseChange.transpose();
 
    //remove parameters which are not given in Hammer
    auto[redBglCentral, redBglCov]  = reduce(bglCentral, bglCov, false);
    //cout << "before diagonalizing" << endl;

    //cout << "Chek reduced correlation matrix!!" << endl;
    //cout << "is symmetric? " << redBglCov.isApprox(redBglCov.transpose(), 1e-8) << endl;;


    // diagonalize cov matrix
    // eigenvectors are columns (see https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/Eigenvalues/EigenSolver.h?ref_type=heads#L171)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(redBglCov);
    Eigen::VectorXd eigenvalues  = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    //cout << " Eigenvector mtarix:" <<endl;
    //printMat(eigenvectors);

    // syntax: ham.setOptions("BtoD*BGL: {avec: [0.00038,0.026905,0.]}");
    paras += "avec: ["  + to_string(redBglCentral(0)) + ", " + to_string(redBglCentral(1)) + ", " + to_string(redBglCentral(2)) + "],";
    paras += "bvec: ["  + to_string(redBglCentral(3)) + ", " + to_string(redBglCentral(4)) + ", " + to_string(redBglCentral(5)) + "],";
    paras += "cvec: ["  + to_string(redBglCentral(6)) + ", " + to_string(redBglCentral(7)) + ", " +                               "],";
    paras += "dvec: ["  + to_string(redBglCentral(8)) + ", " + to_string(redBglCentral(9)) + ", " +                               "],";

    //paras += "avec: ["  +  to_string(0.004605664681641084)    + ", " + to_string(-0.002140593040599278) + ", " + to_string(0.15566982447466055) + "],";
    //paras += "bvec: ["  +  to_string(0.003303529928953319)    + ", " + to_string(-0.004284980385058838) + ", " + to_string(0.17791644334552834) + "],";
    //paras += "cvec: ["  +  to_string(-0.0018867020644757423)  + ", " + to_string(0.022525216948547932)  + ", " +                               "],";
    //paras += "dvec: ["  +  to_string(0.03980443778007538)     + ", " + to_string(-0.1872442367469107)   + ", " +                               "],";

    paras += "}";    

    key = "BstoDs*BGLVar: " + paras;
    //now lets set the variations
    for (size_t i = 1; i < eigenvectors.cols() + 1; ++i){
      //cout << "eigenvalue is: " << eigenvalues(i-1) << ", ";
      for (size_t j = 1; j < eigenvectors.rows() + 1; ++j){
      //cout << "eigenvector entry is " << eigenvectors(j-1,i-1) << endl;
      
        //variations for parameter i in direction up/down * parameter j                                      row j! column i!
        variations["e" + to_string(i)]["up"]["delta_e" + to_string(j)] = sqrt(eigenvalues(i-1)) * eigenvectors(j-1,i-1);
        variations["e" + to_string(i)]["down"]["delta_e" + to_string(j)] = -variations["e" + to_string(i)]["up"]["delta_e" + to_string(j)];
  
      }
    }   


  //cout << "after setting variations" <<endl;

  }
  //set central values via setOptions
  ham.setOptions(key);

  return variations;

}

// get CLN parameters to set via setOptions()
string getCLNParameters(string decay){

  string paras = "{";

  if ((decay == "BsDsMuNu") or (decay == "BsDsTauNu")){

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

  else if ((decay == "BsDs*MuNu") or (decay == "BsDs*TauNu")){

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

const char* getInputFile(const char* signal, const char* file){

  // old
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/29_08_2024_09_46_52/all_signals_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/09_10_2024_14_41_07/dsstartau_flatChunk_0.root"; 
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/08_10_2024_08_30_13/dstau_flatChunk_0.root";
  //const char* fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_10_2024_14_06_09/dsmu_flatChunk_0.root";

  string file_str(file);

  string fin_str;
  const char* fin;

  cout << "reaching this" << endl;
  // use string comparison to compare char* 
  if      (strcmp(signal, "dsmu") == 0) {
    //fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/12_11_2024_18_37_40/*";
    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_06_57_58/" + file_str; //more stats!
    fin = fin_str.c_str();
  }
  else if (strcmp(signal, "dsmu_isgw2") == 0) {
    //fin = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/13_11_2024_07_15_45/*";
    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/17_11_2024_17_26_37/" + file_str; //more stats!
    fin = fin_str.c_str();
  }

  else if (strcmp(signal, "dsstarmu") == 0) {
    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_09_43_21/" + file_str; //more stats!
    fin = fin_str.c_str();
  }

  else if (strcmp(signal, "dstau") == 0) {
    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/15_11_2024_13_50_26/" + file_str; //more stats!
    fin = fin_str.c_str();
  }
  else if (strcmp(signal, "dsstartau") == 0) {
    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/flatNano/16_11_2024_09_45_34/" + file_str; //more stats!
    fin = fin_str.c_str();
  }

  else if (strcmp(signal, "signal") == 0){

    fin_str = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/score_trees/sig_26Sep2024_07h46m21s_cons/" + file_str;
    fin = fin_str.c_str();

  }
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
  cout << "===> Input scheme of " << decay  << " is " << inputScheme << " " << endl;
  ham.addFFScheme("schemeInput", {{hadronicDecay, inputScheme}});
  ham.setFFInputScheme({{hadronicDecay, inputScheme}});

  //define target
  cout << "===> Target scheme of " + decay << " is " << targetScheme << endl;
  ham.addFFScheme( "scheme"+decay , {{hadronicDecay, targetScheme}});

  //////////////////////////////
  // Adapt central CLN values //
  //////////////////////////////
  if ((inputScheme == "CLN") || (targetScheme == "CLN")){

    // if input scheme is CLN, we have to adapt parameters from HQET (EvtGen) -> CLN (Hammer)
    string paras = getCLNParameters(decay);

    // string is of type ham.setOptions("BstoDsBGL: {ChiT: 0.01, ChiL: 0.002}");
    string key;

    if(decay.find("Ds*") != string::npos) key = "BstoDs*CLN: " + paras;
    else                                  key = "BstoDsCLN: "  + paras;
    ham.setOptions(key);

  }

  if ((inputScheme == "BGLVar") || (targetScheme == "BGLVar")){

    
    // if input scheme is BGL, we want to adapt eigenvector naming from a00, .. -> e0, .. 
    setBGLNames(ham, hadronicDecay);   

  }

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
  cout << "target scheme now is " << targetScheme << endl;

  // Set all variations to zero to first calculate the central weights
  for(size_t i = 1; i < pars.size() + 1; i++) { variations_central["delta_e" + to_string(i)] = 0.0;}

  // Set these variations
  ham.setFFEigenvectors(hadronicDecay, targetScheme, variations_central);

  // process event
  ham.processEvent();

  // print central weight
  double weight = ham.getWeight("scheme"+decay);

  std::cout << std::fixed << std::setprecision(10);
  //cout << "central weight " << weight << endl;

  return weight;


};

void getVariationWeights(Hammer::Hammer& ham, vector<string> pars, map<string, map<string, map<string,double>>> vars, string decay, string hadronicDecay, string targetScheme, map<string, double> &w){

  map<string, double> variations;
  vector<string> directions = {"up","down"};

  int counter = 0;
  //Now we loop over every parameter, and calculate its up/down variation weight
  for(size_t i = 1; i < pars.size() + 1; i++){
    for(auto dir: directions){

      //cout << "setting weights for parameter " << name << "in direction " << dir << endl;

      //now for parameter <name> in direction <dir> we set the variations:
      for(size_t j = 1; j < pars.size() + 1; j++){

        variations["delta_e" + to_string(j)] = vars["e" + to_string(i)][dir]["delta_e" + to_string(j)];

      }

      //printVariations(variations);
      ham.setFFEigenvectors(hadronicDecay, targetScheme , variations);

      // process event
      ham.processEvent();

      double weight = ham.getWeight("scheme"+decay);
      //cout << " ====> Gives variation weight delta_" +name+dir + " " << weight <<endl;
 
      //cout << "key is " << "e"+to_string(counter)+"_"+dir << " and set to " << weight << endl; 
      w["e" + to_string(i) + "_"+dir] = weight;

    }

    ++counter;

  }


};


int main(int nargs, char* args[]){
 
  //for printing
  std::cout << std::fixed << std::setprecision(10);

  cout << "Chek correlation matrix!!" << endl;
  cout << "is symmetric? " << latticeCorr.isApprox(latticeCorr.transpose(), 1e-8) << endl;;
  bool isSqCorr = latticeCorr.rows() == latticeCorr.cols();
  cout << "is corr square? " << isSqCorr << endl;
  for (int i = 0; i < latticeCorr.rows(); ++i) {
      for (int j = 0; j < latticeCorr.cols(); ++j) {
          if (std::abs(latticeCorr(i, j) - latticeCorr(j, i)) > 1e-4) {
              std::cout << "Mismatch at (" << i << ", " << j 
                        << "): " << latticeCorr(i, j) << " != " << latticeCorr(j, i) << std::endl;
          }
      }
  }


  cout << "Chek covariance matrix!!" << endl;
  cout << "is Cov symmetric? " << latticeCov.isApprox(latticeCov.transpose(), 1e-8) << endl;;
  bool isSqCov = latticeCov.rows() == latticeCov.cols();
  cout << "is square? " << isSqCov << endl;
  for (int i = 0; i < latticeCov.rows(); ++i) {
      for (int j = 0; j < latticeCov.cols(); ++j) {
          if (std::abs(latticeCov(i, j) - latticeCov(j, i)) > 1e-4) {
              std::cout << "Mismatch at (" << i << ", " << j 
                        << "): " << latticeCov(i, j) << " != " << latticeCov(j, i) << std::endl;
          }
      }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver0(latticeCov);
  Eigen::VectorXd eigenvalues0  = solver0.eigenvalues();
  for(auto eig: eigenvalues0){
    cout << "eigenvalue of initial cov is: " << eig << endl;
  }

  // check if command line arg is given
  if (nargs < 2) {
    //./hammer_all.cc counts already as 1!
    cout << "Please specify signal: dsmu, dstau, dsstarmu, dsstartau";
    exit(1);
  }

  cout << args[0] << endl;
  cout << args[1] << endl;
  cout << args[2] << endl;
  cout << args[3] << endl;
  cout << args[4] << endl;
  cout << args[5] << endl;

  // start setting up hammer
  Hammer::IOBuffer fbBuffer;
  
  // need one hammer instance for every signal!
  Hammer::Hammer hamDsMu;
  Hammer::Hammer hamDsStarMu;
  Hammer::Hammer hamDsTau;
  Hammer::Hammer hamDsStarTau;
 
  //default models 
  string inputDsMu      = "CLN"; 
  string inputDsTau     = "ISGW2"; 
  string inputDsStarMu  = "CLN"; 
  string inputDsStarTau = "ISGW2"; 

  if (string(args[2]) != "default"){
    inputDsMu       = args[2];
    inputDsTau      = args[2];
    inputDsStarMu   = args[2];
    inputDsStarTau  = args[2];

  }

  string target = args[3];
 
  // define input/output FF for every signal
  prepareHammer(hamDsMu,          "BsDsMuNu",       "BsDs",   inputDsMu,       target);
  prepareHammer(hamDsStarMu,      "BsDs*MuNu",      "BsDs*",  inputDsStarMu,   target);
  prepareHammer(hamDsTau,         "BsDsTauNu",      "BsDs",   inputDsTau,      target);
  prepareHammer(hamDsStarTau,     "BsDs*TauNu",     "BsDs*",  inputDsStarTau,  target);
  
  //get parameter and variations
  vector<string> parsDs        = getBGLNames("BsDs");
  vector<string> parsDsStar    = getBGLNames("BsDs*");
  
  string pathDs                = getBGLVariations("BsDs"); 
  string pathDsStar            = getBGLVariations("BsDs*"); 
  
  //import bgl values only for Bs->Ds
  YAML::Node varsDsYaml        = YAML::LoadFile(pathDs);

  //turn into dict for Bs->Ds
  map<string, map<string, map<string,double>>> varsDs;

  //loop over every parameter
  for(size_t i = 1; i < parsDs.size() + 1; i++){
    //loop over up/down weight
    for(size_t j = 1; j < parsDs.size() + 1; j++){
      //write into nested map
      varsDs["e" + to_string(i)]["up"]["delta_e" + to_string(j)]   = varsDsYaml["e" + to_string(i)]["up"]["delta_e" + to_string(j)].as<double>();
      varsDs["e" + to_string(i)]["down"]["delta_e" + to_string(j)] = varsDsYaml["e" + to_string(i)]["down"]["delta_e" + to_string(j)].as<double>();

    }
  }
  
  //load root File for signal given in command line (dsmu, dstau, dsstar mu, dsstar tau)
  TChain* tree = new TChain("tree");
  tree->Add( args[4] ); // give channel and file as argument

  cout << "reached this point savely. " << endl;

  //cout << "tree has entries " << tree->GetEntries() << endl;

  //max nr of events to loop over
  int maxevents = tree->GetEntries();
  //int maxevents = 10;
  if (maxevents > tree->GetEntries()) maxevents = tree->GetEntries();

  //define a place holder for the central_weights
  map<int, double> central_weights;
  //and one for the variational weights (nested map, since we have several parameters)
  map<int, map<string,double>> var_weights;

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

  string dest_str = "/scratch/pahwagne/hammer/" + string(args[1]) + "_" + string(args[3]) + "_" + string(args[5]) + ".root";  
  //string dest_str = string(args[1]) + "_circ.root";  
  const char* dest = dest_str.c_str();

  cout << "destination is " << dest << endl;

  TFile outputFile(dest, "RECREATE");
  TTree* outputTree = tree->CloneTree(0);

  if (access(dest_str.c_str(), W_OK) == 0) {
        std::cout << "Write permission granted for: " << dest_str << std::endl;
    } else {
        std::cerr << "Write permission denied for: " << dest_str << std::endl;
        std::cerr << "Error: " << strerror(errno) << std::endl;
        exit(EXIT_FAILURE); // Exit if permission is not granted
  }

  double central_w;
  double e1_up;
  double e1_down;
  double e2_up;
  double e2_down;
  double e3_up;
  double e3_down;
  double e4_up;
  double e4_down;
  double e5_up;
  double e5_down;
  double e6_up;
  double e6_down;
  double e7_up;
  double e7_down;
  double e8_up;
  double e8_down;
  double e9_up;
  double e9_down;
  double e10_up;
  double e10_down;
  outputTree->Branch("central_w", &central_w, "central_w/D");
  outputTree->Branch("e1_up",     &e1_up,    "e1_up/D");
  outputTree->Branch("e1_down",   &e1_down,  "e1_down/D");
  outputTree->Branch("e2_up",     &e2_up,    "e2_up/D");
  outputTree->Branch("e2_down",   &e2_down,  "e2_down/D");
  outputTree->Branch("e3_up",     &e3_up,    "e3_up/D");
  outputTree->Branch("e3_down",   &e3_down,  "e3_down/D");
  outputTree->Branch("e4_up",     &e4_up,    "e4_up/D");
  outputTree->Branch("e4_down",   &e4_down,  "e4_down/D");
  outputTree->Branch("e5_up",     &e5_up,    "e5_up/D");
  outputTree->Branch("e5_down",   &e5_down,  "e5_down/D");
  outputTree->Branch("e6_up",     &e6_up,    "e6_up/D");
  outputTree->Branch("e6_down",   &e6_down,  "e6_down/D");
  outputTree->Branch("e7_up",     &e7_up,    "e7_up/D");
  outputTree->Branch("e7_down",   &e7_down,  "e7_down/D");
  outputTree->Branch("e8_up",     &e8_up,    "e8_up/D");
  outputTree->Branch("e8_down",   &e8_down,  "e8_down/D");
  outputTree->Branch("e9_up",     &e9_up,    "e9_up/D");
  outputTree->Branch("e9_down",   &e9_down,  "e9_down/D");
  outputTree->Branch("e10_up",    &e10_up,   "e10_up/D");
  outputTree->Branch("e10_down",  &e10_down, "e10_down/D");

  // get max nr of parameters
  int n_vars = max({int(parsDs.size()),int(parsDsStar.size()) });


 
  for(size_t i = 0; i < maxevents; ++i){
  
    if(i%1000 == 0) cout << " ====> Processing event " << i << endl;
  
    //start event
    hamDsMu.initEvent();
    hamDsTau.initEvent();
    hamDsStarMu.initEvent();
    hamDsStarTau.initEvent();
  
    // load entry i into memory 
    tree->GetEntry(i);
  
    //cout << " --------- NEW EVENT of signal " << int(sig) << " -------------" << endl;
    //cout << " bs pt " << bs_pt << " mu pt " << mu_pt << endl;
  
    // Define the decays for which you set up hammer    
    set<int> validDecays = {0,1,10,11};
  
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
    TLorentzVector bsTlv;
    bsTlv.SetPtEtaPhiM(bs_pt,  bs_eta,  bs_phi,  bs_m);
    Hammer::Particle bsPar( Hammer::FourMomentum(bsTlv.E(),  bsTlv.Px(),  bsTlv.Py(),  bsTlv.Pz()),  int(bs_pdgid));
    bs.pt  = bs_pt;     bs.phi  = bs_phi;     bs.eta  = bs_eta;     bs.mass  = bs_m;     bs.pdgid  = bs_pdgid;
  
    //create the process
    Hammer::Process process;

    //placeholders for weights
    double central_weight;
    map<string, double> weights;
    vector<string> directions = {"up","down"}; 
   
    // fill central weight with -1;
    central_w = -1;
    // fill variation weights with -1
    for (size_t i = 0; i < n_vars; i++) {
      for(auto dir:directions){
        weights["e"+to_string(i)+"_"+dir] = -1;
      }
    }

    if (sig == 0){
      //////////////////////////
      // Bs -> Ds + Mu + Nu   //
      //////////////////////////

      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsMu, process, bs, lep, hc);
  
      TLorentzVector dsTlv;
      dsTlv.SetPtEtaPhiM(ds_pt,  ds_eta,  ds_phi,  ds_m);

      if (target == "BGLVar" ){
      
        //variations not used yet 
        auto variations = getBGLSettings(hamDsMu, "BsDsMuNu" , dsTlv, bsTlv);
        //get variation weights
        getVariationWeights(hamDsMu, parsDs, varsDs, "BsDsMuNu", "BstoDs", "BGLVar", weights);

      } 

      // get central weight
      central_weight = getCentralWeight(hamDsMu, parsDs, "BsDsMuNu", "BstoDs", target);
    }
 
    if (sig == 1){
      //////////////////////////
      // Bs -> Ds + Tau+ Nu   //
      //////////////////////////

      lep.pt = tau_pt;     lep.phi = tau_phi;     lep.eta = tau_eta;     lep.mass = tau_m;     lep.pdgid = tau_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsTau, process, bs, lep, hc);
  
      TLorentzVector dsTlv;
      dsTlv.SetPtEtaPhiM(ds_pt,  ds_eta,  ds_phi,  ds_m);
   

      if (target == "BGLVar"){
 
        //variations not used yet 
        auto variations = getBGLSettings(hamDsTau, "BsDsTauNu", dsTlv, bsTlv);
      
        //get variation weights
        getVariationWeights(hamDsTau, parsDs, varsDs, "BsDsTauNu", "BstoDs", "BGLVar", weights);

      } 
      // get central weight
      central_weight = getCentralWeight(hamDsTau, parsDs, "BsDsTauNu", "BstoDs", target);

    }

    if (sig == 10){
      //////////////////////////
      // Bs -> Ds* + Mu+ Nu   //
      //////////////////////////

      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = dsStar_pt;     hc.phi  = dsStar_phi;     hc.eta  = dsStar_eta;     hc.mass  = dsStar_m;     hc.pdgid  = dsStar_pdgid;
  
      //define also DsStar vector for q2 and pPrime calculation
      TLorentzVector dsStarTlv;
      dsStarTlv.SetPtEtaPhiM(dsStar_pt,  dsStar_eta,  dsStar_phi,  dsStar_m);

      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarMu, process, bs, lep, hc);


     //////////////////////////////
     // Adapt central BGL values //
     //////////////////////////////
     if (target == "BGLVar"){
   
       // if input scheme is BGL, we have to adapt central values (done in getBGLSettings)
       // and extract the variations for this event (defined by q2 and pPrime2) and returned
       auto variations = getBGLSettings(hamDsStarMu, "BsDs*MuNu", dsStarTlv, bsTlv);
   
       getVariationWeights(hamDsStarMu, parsDsStar, variations, "BsDs*MuNu", "BstoDs*", "BGLVar", weights);

      } 
 
      // get central weight
      central_weight = getCentralWeight(hamDsStarMu, parsDsStar, "BsDs*MuNu", "BstoDs*", target);

     }

    if (sig == 11){
      //////////////////////////
      // Bs -> Ds* + Tau+ Nu  //
      //////////////////////////

      lep.pt = tau_pt;     lep.phi = tau_phi;     lep.eta = tau_eta;     lep.mass = tau_m;     lep.pdgid = tau_pdgid;
      hc.pt  = dsStar_pt;     hc.phi  = dsStar_phi;     hc.eta  = dsStar_eta;     hc.mass  = dsStar_m;     hc.pdgid  = dsStar_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarTau, process, bs, lep, hc);
 
      //define also DsStar vector for q2 and pPrime calculation
      TLorentzVector dsStarTlv;
      dsStarTlv.SetPtEtaPhiM(dsStar_pt,  dsStar_eta,  dsStar_phi,  dsStar_m);

     if (target == "BGLVar"){
   
       // if input scheme is BGL, we have to adapt central values (done in getBGLSettings)
       // and extract the variations for this event (defined by q2 and pPrime2) and returned
       auto variations = getBGLSettings(hamDsStarTau, "BsDs*TauNu", dsStarTlv, bsTlv);
   
       getVariationWeights(hamDsStarTau, parsDsStar, variations, "BsDs*TauNu", "BstoDs*", "BGLVar", weights);

      } 
 
      // get central weight
      central_weight = getCentralWeight(hamDsStarTau, parsDsStar, "BsDs*TauNu", "BstoDs*", target);

    }



    central_w = central_weight;
    e1_up     = weights["e1_up"]  ; 
    e1_down   = weights["e1_down"]; 
    e2_up     = weights["e2_up"]  ; 
    e2_down   = weights["e2_down"]; 
    e3_up     = weights["e3_up"]  ; 
    e3_down   = weights["e3_down"]; 
    e4_up     = weights["e4_up"]  ; 
    e4_down   = weights["e4_down"]; 
    e5_up     = weights["e5_up"]  ; 
    e5_down   = weights["e5_down"]; 
    e6_up     = weights["e6_up"]  ; 
    e6_down   = weights["e6_down"]; 
    e7_up     = weights["e7_up"]  ; 
    e7_down   = weights["e7_down"]; 
    e8_up     = weights["e8_up"]  ; 
    e8_down   = weights["e8_down"]; 
    e9_up     = weights["e9_up"]  ; 
    e9_down   = weights["e9_down"]; 
    e10_up    = weights["e10_up"]  ; 
    e10_down  = weights["e10_down"]; 

    //cout << "central w is " << central_w << endl;
    //for(int i = 1 ; i < 11 ; ++i){

    //  cout << "variation weight " << i << " is " << weights["e"+to_string(i)+"_up"] << endl;

    //}

    outputTree->Fill(); 
    central_weights[event_int] = central_weight; 
    //var_weights[int(event)]     = weights; 

  } // closing event loop

  outputTree->Write();
  outputFile.Close();

  cout << "saved!" << endl; 
  //cout << "saving to: " << "/scratch/pahwagne/hammer/" + string(args[1]) + "_circ.root" << endl;
  //df.Snapshot("tree", string(args[1]) + "_circ.root");


} //closing main
