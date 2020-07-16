//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
//========================================================= includes
#include "./src/parameters.h"
#include "./src/unpolarised.h"
#include "./src/sivers.h"
#include "./src/transversity.h"
#include "./src/collins.h"

using namespace std;

//========================================================= main
int main(int argc, char **argv)
{

  PARAMETERS P;
  PARTONCONTENT part;

  double Q2 = 2.41;
 
// --- Root
  TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(0);
  plain->SetStatColor(0);
  plain->SetPalette(1);
   

  gROOT->SetStyle("Plain");

  TApplication* theApp = new TApplication("App", &argc, argv);

  DrawSiversDistribution( part,  P, Q2);

  DrawTransversityDistribution( part,  P, Q2);

  DrawCollinsDistribution( part,  P, Q2);

  theApp->Run(true);

  return 0;
}










