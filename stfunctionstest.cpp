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

  cout << " FUTCollins( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  collins << endl;

  cout << "Collins asymmetry = " << collins/unpolarised << endl;

  cout << "Per flavour example:" << endl;

  PARTONCONTENT siversflav = torino10.FUTSiversparton(target, hadron, S, x, z, Q2, PhT);

  cout << "up FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.up << endl;

  cout << "down FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.down << endl;

  cout << "anti_up FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.anti_up << endl;

  cout << "anti_down FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.anti_down << endl;

  cout << "strange FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.strange << endl;

  cout << "anti_strange FUTSivers( " << target << " , " << hadron << " , " <<  S  << " , " <<  x << " , " << z << " , " << Q2 << " , " <<  PhT << ") = " << 
  siversflav.anti_strange << endl;

  return 0;
}