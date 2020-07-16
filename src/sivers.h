//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef __SIVERS_H_
#define __SIVERS_H_
#include "unpolarised.h"
#include "parameters.h"

// --- This is for root use in C++
#include "TApplication.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TPaveLabel.h"
#include "TMultiGraph.h"
#include "TLegend.h" 
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h" 
#include "TF1.h"
#include "TStyle.h" 
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TFile.h" 
#include "TObject.h"
#include "TPaveStats.h"
#include "TText.h"


//**************************************************************************************
//**************************          SIVERS      **************************************
//========================================================= Sivers effect...x dependence
double sivers_x_dependence(double x, double a, double b, double n);

//========================================================= Sivers effect...
void sivers(struct PARTONCONTENT& partcontent, PARAMETERS& Params, double x );


///========================================================= Sivers partcontent for the proton
void SiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double Q2);

///========================================================= Sivers partcontent for all targets
void SiversDistribution(std::string target, struct PARTONCONTENT& partcontent, PARAMETERS& Params, double x, double Q2);

///========================================================= Unpolarised partcontent
void UnpolarisedDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double Q2);

///========================================================= Sivers partcontent KT
void SiversDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double kt, double Q2);

///========================================================= Sivers partcontent First moment
void SiversDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double x, double Q2);

///========================================================= Draw Sivers partcontent
TCanvas* DrawSiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS& Params, double Q2);

#endif //#ifndef __SIVERS_H_