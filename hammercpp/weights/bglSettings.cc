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
#include "bglSettings.h"
#include <Eigen/Dense>

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

tuple<Eigen::MatrixXd, Eigen::MatrixXd> getCorrCovMatrix(bool print){

  int size = 4;

  // Bases is:
  //     V A0 A1 A2
  //  V  x
  // A0  x  x  
  // A1  x  x  x
  // A2  x  x  x  x

  // where the sub-matrix in the lower left corner is 
  // f.e. called a2_v
  // the corresponding vector of variances for variable
  // v is f.e. called v_var

  vector<double> v_var = {
    0.010986180946311711,
    0.27156457203863027,
    0.9982338134983347,
    0.9999954929101635
  }; 

  vector<double> a0_var = {
    0.005666144060237571,
    0.12601100915655922,
    0.9645535514266992,
    0.9989695064536274
  };
  vector<double> a1_var = {
    0.002058935611685443,
    0.053653045313300764,
    0.7745094758048273,
    0.9883994834920944
  };
  vector<double> a2_var = {
    0.011236628965675507,
    0.21619067496595726,
    0.789739961372384,
    0.995181865078575
  };

  vector< vector <double>> v_v = {
    {1.0,       -0.4502,    0.0155,    0.0007512},
    {-0.4502,    1.0,      -0.1043,   -0.003208},
    {0.0155,    -0.1043,    1.0,      -0.0001705},
    {0.0007512, -0.003208, -0.0001705, 1.0}
  };

  vector< vector <double>> a0_v = {
      {-0.002394,  -0.01266,  -0.000435,  -2.199e-05},
      { 0.003344,  0.01268,   0.000509,   2.598e-05},
      {-0.001888,  0.003397,  0.0003099,  1.767e-05},
      { 0.0001465, 0.0001454, 2.376e-05,  1.374e-06}
  };


  vector< vector <double>> a1_v = {
      { 0.008341,  0.009546,  0.0007997,  4.367e-05},
      { 0.01089,   0.02048,   0.0008922,  4.38e-05},
      {-0.000702,  0.01654,   0.0008034,  4.176e-05},
      {-0.001079,  0.0005229, 1.918e-06,  -2.31e-08}
    };


  vector< vector <double>> a2_v = {
      { 0.01155 , -0.025, -0.000755, 3.602e-05},
      { -0.002784, 0.04992, 0.001724 , 8.523e-05},
      {0.003861, -0.002808, 0.0001025, 7.013e-06},
      {0.0008794, 4.2e-05, 2.166e-05, 1.227e-06}
    };

  vector< vector <double>> a0_a0 = {
    {1.0,       -0.3781,    -0.007023,  -0.004275},
    {-0.3781,    1.0,       -0.1681,     0.00345},
    {-0.007023, -0.1681,     1.0,       -0.002318},
    {-0.004275,  0.00345,   -0.002318,   1.0}
  };

  vector< vector <double>> a1_a0 = {
      {0.2776, 0.04791, 0.01723, 0.003371},
      {-0.02225, 0.0006896, 0.03742, -0.008893},
      {0.1044, -0.2142, 0.0506, 0.0279},
      {0.01322, -0.02751, -0.003768, 0.006599}
  };

  vector< vector <double>> a2_a0 = {
      {-0.3816, 0.03317, 0.1246, 0.005469},
      {0.2939, -0.4783, -0.2087, -0.001225},
      {-0.09207, 0.1869, -0.01677, -0.03466},
      {-0.006694, 0.01545, 0.004254, -0.004052}
  };

  vector< vector <double>> a1_a1 = {
    {1.0,       -0.03043,   -0.01583,   -0.01588},
    {-0.03043,   1.0,       -0.3144,    0.02958},
    {-0.01583,  -0.3144,     1.0,       -0.09184},
    {-0.01588,   0.02958,   -0.09184,   1.0} 
  };

  vector< vector <double>> a2_a1 = {
      {0.3327, 0.4326, -0.3016, -0.005564},
      {-0.09657, 0.04303, 0.5087, -0.007639},
      {0.066, -0.1541, 0.485, 0.1149},
      {0.01124, -0.02324, 0.05027, 0.01492}
  };

  vector< vector <double>> a2_a2 = {
    {1.0,       -0.6033,    0.1237,    -0.001915},
    {-0.6033,    1.0,      -0.2094,     0.01039},
    {0.1237,    -0.2094,    1.0,      -0.07082},
    {-0.001915,  0.01039,   -0.07082,   1.0}
  };

  int bigSize = size * size;
  Eigen::MatrixXd corr = Eigen::MatrixXd::Zero(bigSize, bigSize);;
  Eigen::MatrixXd cov  = Eigen::MatrixXd::Zero(bigSize, bigSize);;


  for (size_t i = 0; i < size; ++i){
    for (size_t j = 0; j < size; ++j){

      //diagonal elements 
      corr(i,j)                      = v_v[i][j];
      cov( i,j)                      = corr(i,j) * sqrt(v_var[i]) * sqrt(v_var[j]);

      corr(i + size, j + size)       = a0_a0[i][j];
      cov( i + size, j + size)       = corr(i + size,j + size ) * sqrt(a0_var[i]) * sqrt(a0_var[j]);

      corr(i + 2*size, j + 2*size)   = a1_a1[i][j];
      cov( i + 2*size, j + 2*size)   = corr(i + 2*size,j + 2*size ) * sqrt(a1_var[i]) * sqrt(a1_var[j]);

      corr(i + 3*size, j + 3*size)   = a2_a2[i][j];
      cov( i + 3*size, j + 3*size)   = corr(i + 3*size,j + 3*size ) * sqrt(a2_var[i]) * sqrt(a2_var[j]);

      //off diagonals
      corr(i + size, j)              = a0_v[i][j];
      cov( i + size, j)              = corr(i + size, j) * sqrt(a0_var[i]) * sqrt(v_var[j]);

      corr(i, j + size)              = a0_v[j][i];
      cov( i, j + size)              = corr(i, j + size) * sqrt(a0_var[i]) * sqrt(v_var[j]);

      corr(i + 2*size, j)            = a1_v[i][j];
      cov( i + 2*size, j)            = corr(i + 2*size, j) * sqrt(a1_var[i]) * sqrt(v_var[j]);

      corr(i, j + 2*size)            = a1_v[j][i];
      cov( i, j + 2*size)            = corr(i, j + 2*size) * sqrt(a1_var[i]) * sqrt(v_var[j]);

      corr(i + 3* size, j)           = a2_v[i][j];
      cov( i + 3*size,  j)           = corr(i + 3*size, j) * sqrt(a2_var[i]) * sqrt(v_var[j]);

      corr(i, j + 3* size)           = a2_v[j][i];
      cov( i, j + 3* size)           = corr(i, j + 3*size) * sqrt(a2_var[i]) * sqrt(v_var[j]);

      corr(i + 2*size, j + size)     = a1_a0[i][j];
      cov( i + 2*size, j + size)     = corr(i + 2*size, j + size) * sqrt(a1_var[i]) * sqrt(a0_var[j]);

      corr(i + size, j + 2*size)     = a1_a0[j][i];
      cov( i + size, j + 2*size)     = corr(i + size, j + 2*size) * sqrt(a1_var[i]) * sqrt(a0_var[j]);

      corr(i + 3*size, j + size)     = a2_a0[i][j];
      cov( i + 3*size, j + size)     = corr(i + 3*size, j + size) * sqrt(a2_var[i]) * sqrt(a0_var[j]);

      corr(i + size, j + 3*size)     = a2_a0[j][i];
      cov( i + size, j + 3*size)     = corr(i + size, j + 3*size) * sqrt(a2_var[i]) * sqrt(a0_var[j]);

      corr(i + 3*size, j + 2* size)  = a2_a1[i][j];
      cov( i + 3*size, j + 2* size)  = corr(i + 3*size, j + 2*size) * sqrt(a2_var[i]) * sqrt(a1_var[j]);

      corr(i + 2*size, j + 3* size)  = a2_a1[j][i];
      cov( i + 2*size, j + 3* size)  = corr(i + 2*size, j + 3*size) * sqrt(a2_var[i]) * sqrt(a1_var[j]);

    }
  }

  
  if (print){
    cout << "----- 16x16 CORRELATION MATRIX -----" << endl; 
    for (size_t i = 0; i < corr.rows(); ++i){
      for(size_t j = 0; j < corr.cols(); ++j){
        cout << setw(10) << corr(i,j) << " ";
      }
      cout << endl;
    }
  }
  if (print){
    cout << "----- 16x16 COVARIANCE MATRIX -----" << endl; 
    for (size_t i = 0; i < cov.rows(); ++i){
      for(size_t j = 0; j < cov.cols(); ++j){
        cout << setw(10) << cov(i,j) << " ";
      }
      cout << endl;
    }
  }


  return tuple(corr, cov);

}

Eigen::MatrixXd getBaseChangeMatrix(double q2,double pPrime2, bool print){

  //build matrix which changes basis from V0, A0, A1, A2 to g, f, F1, F2

  int size    = 4;
  int bigSize = 16;

  double mb   = 5.36688; //bs mass
  double md   = 2.1122;  //ds* mass

  // build 10x10 matrix with zeroes
  Eigen::MatrixXd base = Eigen::MatrixXd::Zero(bigSize, bigSize);;

  for (size_t i = 0; i < size; ++i){

    // g
    base(i,i)                  = 2 / (mb + md); 
    // f(
    base(i + size, i + 2*size) = (mb + md);
    // F1(
    base(i + 2*size, i+2*size) = (mb + md) / md * ( - 0.5 * (q2 - mb*mb + md*md ));
    base(i + 2*size, i+3*size) = (mb + md) / md * ( - 2 * mb*mb * pPrime2 /((mb + md)*(mb + md)) );
    // F2(
    base(i + 3*size, i + size) = 2.0;
  }
 
  if (print){
 
    for (size_t i = 0; i < base.rows(); ++i){
      for(size_t j = 0; j < base.cols(); ++j){
        cout << setw(10) << base(i,j) << " ";
      }
      cout << endl;
    }
  }


  return base;
 

}

Eigen::VectorXd getCentralValues(){

  Eigen::VectorXd central(16);
  central << 0.09972089808179423   ,   //V 
             -0.18257567090530047  ,   //V 
             -0.005778560585718485 ,   //V 
             4.1787642909196934e-05,   //V 
                                           
             0.10473977894207792   ,   //A0
             -0.4306269757620783   ,   //A0
             -0.10133044545026366  ,   //A0
             -0.025074536009632522 ,   //A0
                                           
             0.05523807429203482   ,   //A1
             -0.010113200340450334 ,   //A1
             -0.0334447118144628   ,   //A1
             0.06348503204545279   ,   //A1
                                           
             0.05909940832954398   ,   //A2
             -0.1087887618586812   ,   //A2
             -0.24726887556559585  ,   //A2
             -0.04689433559347342;     //A2

  return central;                     

}


Eigen::MatrixXd getCovariance(Eigen::MatrixXd corr, Eigen::VectorXd sig){


  int bigSize = corr.rows();
  Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(bigSize, bigSize); 

  for (size_t i = 0; i < corr.rows(); ++i){
    for(size_t j = 0; j < corr.cols(); ++j){

      cov(i,j) = corr(i,j) * sig(i) * sig(j);

    }
  }
  
  return cov;

}

tuple<Eigen::VectorXd, Eigen::MatrixXd> reduce(Eigen::VectorXd vec, Eigen::MatrixXd mat, bool print){

  if (print){
 
    cout << "----- 16x16 MATRIX -----" << endl;
    for (size_t i = 0; i < mat.rows(); ++i){
      for(size_t j = 0; j < mat.cols(); ++j){
        cout << setw(14) << mat(i,j) << " ";
      }
      cout << endl;
    }
  }

  if (print){
 
    cout << "----- 16 VECTOR -----" << endl;
    for (size_t i = 0; i < mat.rows(); ++i){
      cout << setw(14) << vec(i) << endl;
    }
  }

  int bigSize = 16;
  vector<int> skip = {3,7,8,11,14,15};

  // now remove these cols and rows
  Eigen::MatrixXd redVec = Eigen::VectorXd::Zero(10);
  Eigen::MatrixXd redMat = Eigen::MatrixXd::Zero(10, 10);

  int new_row = 0;
  for (size_t i = 0; i < bigSize; ++i){
    
    int new_col = 0;

    //skip rows  in skip
    if(find(skip.begin(), skip.end(), i) != skip.end()) continue;

    redVec(new_row) = vec(i);

    for (size_t j = 0; j < bigSize; ++j){

      //skip cols in skip
      if(find(skip.begin(), skip.end(), j) != skip.end()) continue;
  
      redMat(new_row, new_col) = mat(i,j); 
      ++new_col;
    }
    ++new_row;
  }

  if (print){
 
    cout << "----- 10x10 REDUCED MATRIX -----" << endl;
    for (size_t i = 0; i < redMat.rows(); ++i){
      for(size_t j = 0; j < redMat.cols(); ++j){
        cout << setw(14) << redMat(i,j) << " ";
      }
      cout << endl;
    }
  }

  if (print){
 
    cout << "----- 10 VECTOR -----" << endl;
    for (size_t i = 0; i < redMat.rows(); ++i){
        cout << setw(14) << redVec(i) << endl;
    }
  }

 return tuple(redVec, redMat);

}

void printMat(Eigen::MatrixXd mat){

    for (size_t i = 0; i < mat.rows(); ++i){
      for(size_t j = 0; j < mat.cols(); ++j){
        cout << setw(14) << mat(i,j) << endl;
      }
    }
}



int main() {

    // get  lattice central values and corr/cov matrix
    Eigen::VectorXd latticeCentral = getCentralValues();
    auto[latticeCorr, latticeCov]  = getCorrCovMatrix(true);

    // get change of basis
    Eigen::MatrixXd baseChange     = getBaseChangeMatrix(4, 2, true);

    //transform the central values!
    Eigen::VectorXd bglCentral  = baseChange * latticeCentral; 
    //transform covariance matrix - it is a bilinear form!
    Eigen::MatrixXd bglCov      = baseChange * latticeCov * baseChange.transpose(); 

    // remove parameters which are not given in Hammer
    auto[redBglCentral, redBglCov]  = reduce(bglCentral, bglCov, true); 

    // diagonalize cov matrix
    // eigenvectors are columns (see https://gitlab.com/libeigen/eigen/-/blob/master/Eigen/src/Eigenvalues/EigenSolver.h?ref_type=heads#L171)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(redBglCov); 
    Eigen::VectorXd eigenvalues  = solver.eigenvalues(); 
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();  

    cout << "Eigenvalues:\n" << eigenvalues << std::endl;
    cout << "Eigenvectors (columns):\n" << eigenvectors << std::endl; 

    return 0;

}

