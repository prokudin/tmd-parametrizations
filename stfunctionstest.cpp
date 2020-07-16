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
  double x  = 0.2;
  double z  = 0.4;
  double Q2 = 2.41;
  double PhT = 0.6;

  double unpolarised = torino10.FUU(target, hadron, S, x, z, Q2, PhT);
  double sivers = torino10.FUTSivers(target, hadron, S, x, z, Q2, PhT);
  double collins = torino10.FUTCollins(target, hadron, S, x, z, Q2, PhT);
 
  cout << " FUU( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  unpolarised << endl;

  cout << " FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  sivers << endl;

  cout << "Sivers asymmetry = " << sivers/unpolarised << endl;

  cout << " FUTCollins( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  collins << endl;

  cout << "Collins asymmetry = " << collins/unpolarised << endl;

  return 0;
}