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
vector<string> getBGLNames(string channel){

  vector<string> pars;

  if (channel == "BsDs") pars = {"a00","a01","a02","a03","ap0","ap1","ap2","ap3"}; //scalar
  else                   pars = {"a0","a1","a2","b0","b1","b2","c1","c2", "d0", "d1"}; //vectorial

  return pars; 
}

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

// get BGL variations
string getBGLVariations(string channel){

  string path;
  if (channel == "BsDs") path = "/work/pahwagne/RDsTools/hammer/bgl_scalar_variations.yaml";
  else                   path = "/work/pahwagne/RDsTools/hammer/bgl_vector_variations.yaml";

  return path;
}

// get BGL parameters to set via setOptions()

string getBGLParameters(string decay){

  string paras = "{";

  if ((decay == "BsDsMuNu") or (decay == "BsDsTauNu")){
    cout << "Setting scalar BGL central values" << endl;

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
    //paras += "}";    

    paras += "a0: ["  +  to_string(0.052255946347001495)    + ", " + to_string(-0.16027634967890908) + ", " + to_string(0.014141836205563255) + ", " + to_string(0.0) + "],";
    paras += "ap: ["  +  to_string(0.0017893827864468802)   + ", " + to_string(-0.004691380424494185) + ", " + to_string( -0.015708534616906505) + ", " + to_string(0.0) + "],";
    paras += "Chim: 0.007168,";
    paras += "Chip: 0.012539,";
    paras += "ChimL: 0.025042,";
    paras += "}";    


  }

  else{
    cout << "Setting vector BGL central values" << endl;

    // TODO: import from yaml file
    double a0 = 0.004605664681641084; 
    double a1 = -0.002140593040599278; 
    double a2 = 0.15566982447466055;
    double b0 = 0.003303529928953319; 
    double b1 = -0.004284980385058838;
    double b2 = 0.17791644334552834;
    //double c0 = missing in hammer
    double c1 = -0.0018867020644757423;  
    double c2 = 0.022525216948547932;
    double d0 = 0.03980443778007538;
    double d1 = -0.1872442367469107; 
    //double d2 = missing in hammer

    //paras += "a0: " + to_string(a0) + ", ";
    //paras += "a1: " + to_string(a1) + ", ";
    //paras += "a2: " + to_string(a2) + ", ";
    //paras += "b0: " + to_string(b0) + ", ";
    //paras += "b1: " + to_string(b1) + ", ";
    //paras += "b2: " + to_string(b2) + ", ";
    //paras += "c1: " + to_string(c1) + ", ";
    //paras += "c2: " + to_string(c2) + ", ";
    //paras += "d0: " + to_string(d0) + ", ";
    //paras += "d1: " + to_string(d1) + ", ";

    //paras += "a0: " + to_string(0.0) + ", ";
    //paras += "a1: " + to_string(0.0) + ", ";
    //paras += "a2: " + to_string(0.0) + ", ";
    //paras += "b0: " + to_string(0.0) + ", ";
    //paras += "b1: " + to_string(0.0) + ", ";
    //paras += "b2: " + to_string(0.0) + ", ";
    //paras += "c1: " + to_string(0.0) + ", ";
    //paras += "c2: " + to_string(0.0) + ", ";
    //paras += "d0: " + to_string(0.0) + ", ";
    //paras += "d1: " + to_string(0.0) + ", ";


    //paras += "}";    

    // Cohen paper
    //paras += "avec: ["  +  to_string(0.004605664681641084)    + ", " + to_string(-0.002140593040599278) + ", " + to_string(0.15566982447466055) + "],";
    //paras += "bvec: ["  +  to_string(0.003303529928953319)    + ", " + to_string(-0.004284980385058838) + ", " + to_string(0.17791644334552834) + "],";
    //paras += "cvec: ["  +  to_string(-0.0018867020644757423)  + ", " + to_string(0.022525216948547932)  + ", " +                               "],";
    //paras += "dvec: ["  +  to_string(0.03980443778007538)     + ", " + to_string(-0.1872442367469107)   + ", " +                               "],";

    //paras += "Chim: 0.007168,";
    //paras += "Chip: 0.012539,";
    //paras += "ChimL: 0.025042,";

    //paras += "Vcb: 0.99344,";
   
    //Olmo
    //paras += "avec: ["  +  to_string(0.03200)    + ", " + to_string(-0.14800)   + ", " + to_string(-0.60000) + "],";
    //paras += "bvec: ["  +  to_string(0.01246)    + ", " + to_string(0.00380)    + ", " + to_string(0.02000)  + "],";
    //paras += "cvec: ["  +  to_string(0.00008)    + ", " + to_string(0.08000)    + ", " +                       "],";
    //paras += "dvec: ["  +  to_string(0.05260)    + ", " + to_string(-0.19400)   + ", " +                       "],";
 
    //taken from table XXXXII in paper
    //paras += "avec: [0.0300, 0.0, -0.04],"; 
    //paras += "bvec: [0.01420, -0.019, -0.0],";
    //paras += "cvec: [-0.0018, -0.041],";
    //paras += "dvec: [0.0384, -0.077]";

    //converted from paper
    //paras += "avec: [0.026667, -0.048823, -0.001545],"; 
    //paras += "bvec: [0.413130, -0.075637, -0.250136],";
    //paras += "cvec: [1.206462, 2.327528],";
    //paras += "dvec: [0.209480, -0.861254]";

    //Judd Harrison 
    //paras += "avec: [0.0280,  -0.056, 0.08],";    //g
    //paras += "bvec: [0.01257, 0.010, -0.32],";    //f
    //paras += "cvec: [-0.0042, -0.066],";          //F1
    //paras += "dvec: [0.0487, -0.240],";           //F2
    //paras += "Vcb:  0.99344,";
    //paras += "}";    

    paras += "avec: [0.0282,  -0.072, 0.5],";    //g
    paras += "bvec: [0.01256, 0.012, -0.40],";    //f
    paras += "cvec: [-0.0038, -0.090],";          //F1
    paras += "dvec: [0.0480, -0.17],";           //F2
    paras += "Vcb:  0.99344,";
    paras += "}";    


  }

  return paras;

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

  //////////////////////////////
  // Adapt central BGL values //
  //////////////////////////////
  if ((inputScheme == "BGLVar") || (targetScheme == "BGLVar")){

    // if input scheme is CLN, we have to adapt parameters from HQET (EvtGen) -> CLN (Hammer)
    string paras = getBGLParameters(decay);

    // string is of type ham.setOptions("BstoDsBGL: {ChiT: 0.01, ChiL: 0.002}");
    string key;

    if(decay.find("Ds*") != string::npos) key = "BstoDs*BGLVar: " + paras;
    else                                  key = "BstoDsBGLVar: "  + paras;
    ham.setOptions(key);
 

    //// ################# BGL for D* ##############################
    //// Other options still from https://arxiv.org/pdf/2105.14019.pdf
    //string centralValuesOpt = "BstoDs*BGLVar: {";
    ////  Tab 9 line 2: f,F1 (1+)
    //centralValuesOpt += Form("BcStatesf: [6.739, 6.750, 7.145, 7.150], ");
    //// Hammer default                     6.730, 6.736, 7.135, 7.142}; //GeV
    ////  Tab 9 line 1: g (1-)
    //centralValuesOpt += Form("BcStatesg: [6.329, 6.920, 7.020], ");
    //// Hammer default                     6.337, 6.899, 7.012, 7.280}; //GeV
    ////  Tab 9 line 3: F2 (0-)
    //centralValuesOpt += Form("BcStatesP1: [6.275, 6.842, 7.250], ");
    //// Hammer default                      6.275, 6.842, 7.250}; //GeV
    //// Vcb from abstract rest from Tab 10 (chi+ and - inverted)
    //// centralValuesOpt += Form("Vcb: 38.4e-3, Chim: 3.894e-4, Chip: 5.131e-4, ChimL: 19.42e-3");
    //// Setting Vcb at 1/eta_EW = 1/1.0066. = 0.99344
    ////centralValuesOpt += Form("Vcb: 0.99344, Chim: 3.894e-4, Chip: 5.131e-4, ChimL: 19.42e-3");
    //centralValuesOpt += Form("Chim: 3.894e-4, Chip: 5.131e-4, ChimL: 19.42e-3");
    //// Hammer defaults        Vcb: 41.5e-3, Chim: 3.068e-4, Chip: 5.280e-4, ChimL:  2.466e-3
    //centralValuesOpt += "}";
    //ham.setOptions(centralValuesOpt);

    //only fpr judd! 
    //double chim = 0.00737058/(4.78*4.78);
    //double chip = 0.0126835/(4.78*4.78);
    //double chiml = 0.0249307;
    //double mb = 4.78;
    //double mc = mb - 3.4;
    //string outer = "";
    //outer += "{Chim: " + to_string(chim) +","; 
    //outer += "Chip: " + to_string(chip) +","; 
    //outer += "ChiL: " + to_string(chiml) +",";
    //outer += "BcStatesf:  [6.745, 6.75, 7.15, 7.15], ";
    //outer += "BcStatesg:  [6.335, 6.926, 7.02 , 7.28], ";
    //outer += "BcStatesP1: [6.275, 6.872, 7.25], "; 
    //outer += "mb: " + to_string(mb) + ", ";
    //outer += "mc: " + to_string(mc) + ", ";
    //outer += "Vcb: 0.99344";
    //outer += "}";
    //cout << outer << endl;
    //ham.setOptions("BstoDs*BGLVar: " + outer);


    // if input scheme is BGL, we want to adapt parameter naming from a00, .. -> e0, .. 
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
  //cout << "lepton " << l.pt << l.eta << l.phi << endl;
  //cout << "b " << b.pt << b.eta << b.phi << endl;
  //cout << "h " << h.pt << h.eta << h.phi << endl;

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

  //////cout << "target scheme now is " << targetScheme << endl;
  map<string, double> variations_central;

  // Set all variations to zero to first calculate the central weights
  for(size_t i = 1; i < pars.size() + 1; i++) { variations_central["delta_e" + to_string(i)] = 0.0;}

  // Set these variations
  ham.setFFEigenvectors(hadronicDecay, targetScheme, variations_central);

  // process event
  ham.processEvent();

  // print central weight
  double weight = ham.getWeight("scheme"+decay);
  float weightf = ham.getWeight("scheme"+decay);

  //cout << "central weight double " << weight << endl;
  //cout << "central weight float" << weightf << endl;

  return weight;


};

void getVariationWeights(Hammer::Hammer& ham, vector<string> pars, YAML::Node vars, string decay, string hadronicDecay, string targetScheme, map<string, double> &w){

  map<string, double> variations;
  vector<string> directions = {"up","down"};

  int counter = 0;
  //Now we loop over every parameter, and calculate its up/down variation weight
  for(size_t i = 1; i < pars.size() + 1; i++){
    for(auto dir: directions){

      //cout << "setting weights for parameter " << name << "in direction " << dir << endl;

      //now for parameter <name> in direction <dir> we set the variations:
      for(size_t j = 1; j < pars.size() + 1; j++){

        variations["delta_e" + to_string(j)] = vars["e" + to_string(i)][dir]["delta_e" + to_string(j)].as<double>();

      }

      //printVariations(variations);

      ham.setFFEigenvectors(hadronicDecay, targetScheme, variations);

      // process event
      ham.processEvent();

      double weight = ham.getWeight("scheme"+decay);
      //cout << " ====> Gives variation weight delta_" +name+dir + " " << weight <<endl;
 
      //cout << "key is " << "e"+to_string(counter)+"_"+dir << " and set to " << weight << endl; 
      w["e"+to_string(i)+"_"+dir] = weight;

    }

    ++counter;

  }


};

int main(int nargs, char* args[]){

  //for printing
  std::cout << std::fixed << std::setprecision(10);

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
 
  //
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
  prepareHammer(hamDsMu,          "BsDsMuNu",       "BsDs",   inputDsMu,       args[3]);
  prepareHammer(hamDsStarMu,      "BsDs*MuNu",      "BsDs*",  inputDsStarMu,   args[3]);
  prepareHammer(hamDsTau,         "BsDsTauNu",      "BsDs",   inputDsTau,      args[3]);
  prepareHammer(hamDsStarTau,     "BsDs*TauNu",     "BsDs*",  inputDsStarTau,  args[3]);
  
  //get parameter and variations
  vector<string> parsDs        = getBGLNames("BsDs");
  vector<string> parsDsStar    = getBGLNames("BsDs*");
  
  string pathDs                = getBGLVariations("BsDs"); 
  string pathDsStar            = getBGLVariations("BsDs*"); 
  
  YAML::Node varsDs            = YAML::LoadFile(pathDs);
  YAML::Node varsDsStar        = YAML::LoadFile(pathDsStar);
  
  //load root File for signal given in command line (dsmu, dstau, dsstar mu, dsstar tau)
  TChain* tree = new TChain("tree");
  tree->Add( args[4] ); // give channel and file as argument

  cout << "reached this point savely. " << endl;

  //cout << "tree has entries " << tree->GetEntries() << endl;

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


  float central_w;
  float e1_up;
  float e1_down;
  float e2_up;
  float e2_down;
  float e3_up;
  float e3_down;
  float e4_up;
  float e4_down;
  float e5_up;
  float e5_down;
  float e6_up;
  float e6_down;
  float e7_up;
  float e7_down;
  float e8_up;
  float e8_down;
  float e9_up;
  float e9_down;
  float e10_up;
  float e10_down;
  outputTree->Branch("central_w", &central_w, "central_w/F");
  outputTree->Branch("e1_up",   &e1_up,   "e1_up/F");
  outputTree->Branch("e1_down", &e1_down, "e1_down/F");
  outputTree->Branch("e2_up",   &e2_up,   "e2_up/F");
  outputTree->Branch("e2_down", &e2_down, "e2_down/F");
  outputTree->Branch("e3_up",   &e3_up,   "e3_up/F");
  outputTree->Branch("e3_down", &e3_down, "e3_down/F");
  outputTree->Branch("e4_up",   &e4_up,   "e4_up/F");
  outputTree->Branch("e4_down", &e4_down, "e4_down/F");
  outputTree->Branch("e5_up",   &e5_up,   "e5_up/F");
  outputTree->Branch("e5_down", &e5_down, "e5_down/F");
  outputTree->Branch("e6_up",   &e6_up,   "e6_up/F");
  outputTree->Branch("e6_down", &e6_down, "e6_down/F");
  outputTree->Branch("e7_up",   &e7_up,   "e7_up/F");
  outputTree->Branch("e7_down", &e7_down, "e7_down/F");
  outputTree->Branch("e8_up",   &e8_up,   "e8_up/F");
  outputTree->Branch("e8_down", &e8_down, "e8_down/F");
  outputTree->Branch("e9_up",   &e9_up,   "e9_up/F");
  outputTree->Branch("e9_down", &e9_down, "e9_down/F");
  outputTree->Branch("e10_up",   &e10_up,   "e10_up/F");
  outputTree->Branch("e10_down", &e10_down, "e10_down/F");

  // get max nr of parameters
  int n_vars = max({int(parsDs.size()),int(parsDsStar.size()) });
 
  for(size_t i = 0; i < maxevents; ++i){
  
    if(i%100 == 0) cout << " ====> Processing event " << i << endl;
  
    //start event
    hamDsMu.initEvent();
    hamDsStarMu.initEvent();
    hamDsTau.initEvent();
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
    ROOT::Math::PtEtaPhiMVector bsTlv (bs_pt,  bs_eta,  bs_phi,  bs_m);
    Hammer::Particle bsPar( Hammer::FourMomentum(bsTlv.e(),  bsTlv.px(),  bsTlv.py(),  bsTlv.pz()),  int(bs_pdgid));
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
  
      // get central weight
      central_weight = getCentralWeight(hamDsMu, parsDs, "BsDsMuNu", "BstoDs", args[3]);

      if (strcmp(args[3],"BGLVar") == 0){
       
        //get variation weights
        getVariationWeights(hamDsMu, parsDs, varsDs, "BsDsMuNu", "BstoDs", "BGLVar", weights);

      } 
 
      //cout << "At index: " << i << " ====> saving " << central_weight << " for event nr " << event_int << endl; 
      //cout << "At index: " << i << " ====> saving " << weights["e2_up"] << " variation for event nr " << event_int << endl; 
    }
 
    if (sig == 1){
      //////////////////////////
      // Bs -> Ds + Tau+ Nu   //
      //////////////////////////

      lep.pt = tau_pt;     lep.phi = tau_phi;     lep.eta = tau_eta;     lep.mass = tau_m;     lep.pdgid = tau_pdgid;
      hc.pt  = ds_pt;     hc.phi  = ds_phi;     hc.eta  = ds_eta;     hc.mass  = ds_m;     hc.pdgid  = ds_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsTau, process, bs, lep, hc);
  
      // get central weight
      central_weight = getCentralWeight(hamDsTau, parsDs, "BsDsTauNu", "BstoDs", args[3]);

      if (strcmp(args[3],"BGLVar") == 0){
       
        //get variation weights
        getVariationWeights(hamDsTau, parsDs, varsDs, "BsDsTauNu", "BstoDs", "BGLVar", weights);

      } 
 
      //cout << "At index: " << i << " ====> saving " << central_weight << " for event nr " << event_int << endl; 
      //cout << "At index: " << i << " ====> saving " << weights["e2_up"] << " variation for event nr " << event_int << endl; 
    }

    if (sig == 10){
      //////////////////////////
      // Bs -> Ds* + Mu+ Nu   //
      //////////////////////////

      lep.pt = mu_pt;     lep.phi = mu_phi;     lep.eta = mu_eta;     lep.mass = mu_m;     lep.pdgid = mu_pdgid;
      hc.pt  = dsStar_pt;     hc.phi  = dsStar_phi;     hc.eta  = dsStar_eta;     hc.mass  = dsStar_m;     hc.pdgid  = dsStar_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarMu, process, bs, lep, hc);
  
      if (strcmp(args[3],"BGLVar") == 0){
       
        //get variation weights
        getVariationWeights(hamDsStarMu, parsDsStar, varsDsStar, "BsDs*MuNu", "BstoDs*", "BGLVar", weights);

      } 
      // get central weight
      central_weight = getCentralWeight(hamDsStarMu, parsDsStar, "BsDs*MuNu", "BstoDs*", args[3]);

      //cout << "At index: " << i << " ====> saving " << central_weight << " for event nr " << event_int << endl; 
      //cout << "At index: " << i << " ====> saving " << weights["e2_up"] << " variation for event nr " << event_int << endl; 
    
     }

    if (sig == 11){
      //////////////////////////
      // Bs -> Ds* + Tau+ Nu  //
      //////////////////////////

      lep.pt = tau_pt;     lep.phi = tau_phi;     lep.eta = tau_eta;     lep.mass = tau_m;     lep.pdgid = tau_pdgid;
      hc.pt  = dsStar_pt;     hc.phi  = dsStar_phi;     hc.eta  = dsStar_eta;     hc.mass  = dsStar_m;     hc.pdgid  = dsStar_pdgid;
  
      // define decay chain, particle p4 and start process 
      defineDecay(hamDsStarTau, process, bs, lep, hc);
  
      // get central weight
      central_weight = getCentralWeight(hamDsStarTau, parsDsStar, "BsDs*TauNu", "BstoDs*", args[3]);

      if (strcmp(args[3],"BGLVar") == 0){
       
        //get variation weights
        getVariationWeights(hamDsStarTau, parsDsStar, varsDsStar, "BsDs*TauNu", "BstoDs*", "BGLVar", weights);

      } 
 
      //cout << "At index: " << i << " ====> saving " << central_weight << " for event nr " << event_int << endl; 
      //cout << "At index: " << i << " ====> saving " << weights["e2_up"] << " variation for event nr " << event_int << endl; 
    }

    //cout << "central w is " << central_weight << endl;

    //for(int i = 1 ; i < 11 ; ++i){

    //  cout << "variation weight " << i << " is " << weights["e"+to_string(i)+"_up"] << endl;

    //}


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


    outputTree->Fill(); 
    central_weights[event_int] = central_weight; 
    //var_weights[int(event)]     = weights; 

  } // closing event loop

  outputTree->Write();
  outputFile.Close();
 
  //cout << "saving to: " << "/scratch/pahwagne/hammer/" + string(args[1]) + "_circ.root" << endl;
  //df.Snapshot("tree", string(args[1]) + "_circ.root");


} //closing main

