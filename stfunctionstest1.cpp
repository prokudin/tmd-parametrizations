//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
//========================================================= includes
#include "./src/parameters.h"
// #include "./src/unpolarised.h"
// #include "./src/sivers.h"
// #include "./src/transversity.h"
#include "./src/stfunctions.h"
#include <iostream>

using namespace std;

//========================================================= main
int main(int argc, char **argv)
{
  TMD torino10;
  string target = "proton";
  string hadron = "pi+";

  double S  = 17.*17.;
  double x  = 0.12;
  double z  = 0.32;
  double Q2 = 2.41;
  double PhT = 0.14;

  double unpolarised = torino10.FUU(target, hadron, S, x, z, Q2, PhT);
  double sivers = torino10.FUTSivers(target, hadron, S, x, z, Q2, PhT);
//  double collins = torino10.FUTCollins(target, hadron, S, x, z, Q2, PhT);

//    FUU( proton , pi+ , 289 , 0.2 , 0.4 , 2.41 , 0.6) = 1.25882
//  FUTSivers( proton , pi+ , 289 , 0.2 , 0.4 , 2.41 , 0.6) = 0.111049
// Sivers asymmetry = 0.088217
//  FUTCollins( proton , pi+ , 289 , 0.2 , 0.4 , 2.41 , 0.6) = 0.0058584
// Collins asymmetry = 0.00465389
 
  cout << " FUU( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  unpolarised << endl;

  cout << " FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  sivers << endl;

  cout << "Sivers asymmetry = " << sivers/unpolarised << endl;

  int N = 100;
  double Q2min = 1.5;
  double Q2max = 130;
  double dQ2 = (Q2max - Q2min)/N;
    
  for (int i = 0; i<N; i++){
      Q2 = Q2min + i*dQ2;
      unpolarised = torino10.FUU(target, hadron, S, x, z, Q2, PhT);
      sivers = torino10.FUTSivers(target, hadron, S, x, z, Q2, PhT);
      cout << Q2 << " "<< sivers/unpolarised << endl;
  }

  return 0;
}