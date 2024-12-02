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
#include "calculateCovMatrix.h"

////////////////////////////////////////////////////////////////////////////////////////////
// Command line args:                                                                     //
// args[1] = dsmu, ... (i.e. the signal we want to calculate the corr matrix for )        //
////////////////////////////////////////////////////////////////////////////////////////////


using namespace std;

void boostDsStar(TLorentzVector&  dsStar, TLorentzVector bs){

  //boost dsStar into Bs rest frame
  TVector3 bBoost = bs.BoostVector();
  dsStar.Boost(-bBoost);

}

void buildCorrMatrix(){

  int size = 4;

  // Bases is:
  //     V A0 A1 A2
  //  V  x
  // A0  x  x  
  // A1  x  x  x
  // A2  x  x  x  x

  // where the sub-matrix in the lower left corner is 
  // f.e. called a2_v

  double v_v[size][size] = {
      {1.0,       0.0,        0.0 ,        0.0},
      {-0.4502,   1.0,        0.0,         0.0 },
      {0.0155,    -0.1043,    1.0,         0.0},
      {0.0007512, -0.003208,  -0.0001705,  1.0}
  };

  double a0_v[size][size] = {
      {-0.002394,  -0.01266,  -0.000435,  -2.199e-05},
      { 0.003344,  0.01268,   0.000509,   2.598e-05},
      {-0.001888,  0.003397,  0.0003099,  1.767e-05},
      { 0.0001465, 0.0001454, 2.376e-05,  1.374e-06}
  };


  double a1_v[size][size]{
      { 0.008341,  0.009546,  0.0007997,  4.367e-05},
      { 0.01089,   0.02048,   0.0008922,  4.38e-05},
      {-0.000702,  0.01654,   0.0008034,  4.176e-05},
      {-0.001079,  0.0005229, 1.918e-06,  -2.31e-08}
    };


  double a2_v[size][size]{
      { 0.01155 , -0.025, -0.000755, 3.602e-05},
      { -0.002784, -0.002784, 0.001724 , 8.523e-05},
      {0.003861, -0.002808, 0.0001025, 7.013e-06},
      {0.0008794, 4.2e-05, 2.166e-05, 1.227e-06}
    };

  double a0_a0[size][size]{

      {1.0,     0.0,    0.0,    0.0,      },
      {-0.3781, 1.0,    0.0,    0.0       },
      {-0.007023, -0.1681, 1.0, 0.0       },
      {-0.004275, 0.00345, -0.002318, 1.0,}
 
  };
  double a1_a0[size][size]{

      {0.2776, 0.04791, 0.01723, 0.003371},
      {-0.02225, 0.0006896, 0.03742, -0.008893},
      {0.1044, -0.2142, 0.0506, 0.0279},
      {0.01322, -0.02751, -0.003768, 0.006599}

  };

  double a2_a0[size][size]{

      {-0.3816, 0.03317, 0.1246, 0.005469},
      {0.2939, -0.4783, -0.2087, -0.001225},
      {-0.09207, 0.1869, -0.01677, -0.03466},
      {-0.006694, 0.01545, 0.004254, -0.004052}
 
  };
  double a1_a1[size][size]{

      {1.0,    , 0.0,     0.0,       0.0},
      {-0.03043, 1.0,     0.0,       0.0},
      {-0.01583, -0.3144,  1.0,      0.0},
      {-0.01588, 0.02958,  -0.09184, 1.0}
 
  };

  double a2_a1[size][size]{

      {0.3327, 0.4326, -0.3016, -0.005564},
      {-0.09657, 0.04303, 0.5087, -0.007639},
      {0.066, -0.1541, 0.485, 0.1149},
      {0.01124, -0.02324, 0.05027, 0.01492}
  };

  double a2_a2[size][size]{

      {1.0,     0.0,     0.0,    0.0 },
      {-0.6033, 1.0,     0.0,    0.0 },
      {0.1237, -0.2094, 1.0,     0.0},
      {-0.001915, 0.01039, -0.07082, 1.0}
 
  };

  bigSize = size * size;
  double corrMatrix[bigSize][bigSize] = {0};

  for (size_t i; i < size; ++i){
    for (size_t j; j <= i; ++j){

      //diagonal elements 
      corrMatrix[i][j]                     = v_v[i][j];
      corrMatrix[j][i]                     = v_v[j][i];

      corrMatrix[i + size][j + size]       = a0_a0[i][j];
      corrMatrix[j + size][i + size]       = a0_a0[j][i];

      corrMatrix[i + 2*size][j + 2*size]   = a1_a1[i][j];
      corrMatrix[j + 2*size][i + 2*size]   = a1_a1[j][i];

      corrMatrix[i + 3*size][j + 3*size]   = a2_a2[i][j];
      corrMatrix[j + 3*size][i + 3*size]   = a2_a2[j][i];

      //off diagonals
      corrMatrix[i][j + size]              = a0_v[i][j];
      corrMatrix[j][i + size]              = a0_v[j][i];

      corrMatrix[i][j + 2* size]           = a1_v[i][j];
      corrMatrix[j][i + 2* size]           = a1_v[j][i];

      corrMatrix[i][j + 3* size]           = a2_v[i][j];
      corrMatrix[j][i + 3* size]           = a2_v[j][i];

      //off diagonals
      corrMatrix[i + size][j + 2* size]    = a1_a0[i][j];
      corrMatrix[j + size][i+ 2* size]    = a1_a0[j][i];

      corrMatrix[i + size][j + 3* size]    = a2_a0[i][j];
      corrMatrix[j + size][i + 3* size]    = a2_a0[j][i];

      //off diagonals
      corrMatrix[i + 2*size][j + 3* size]  = a2_a1[i][j];
      corrMatrix[j + 2*size][i + 3* size]  = a2_a1[j][i];


  }

}



