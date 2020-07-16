//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef __COLLINS_H_
#define __COLLINS_H_
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
//**************************          COLLINS     **************************************
//========================================================= Collins effect...z dependence
double collins_z_dependence(double z, double a, double b, double n);

//========================================================= Fragmentation PI^+
void Fragmentation(struct PARTONCONTENT& fragmentation, double z, double Q2);

//========================================================= Collins effect PI^+...
void collins(struct PARTONCONTENT& fragmentation, PARAMETERS& Params, double z);

//========================================================= Collins partcontent
void CollinsDistribution( struct PARTONCONTENT& fragmentation, PARAMETERS& Params, double z, double Q2);

void CollinsFragmentation( std::string hadron, struct PARTONCONTENT& fragment, PARAMETERS& Params, double z, double Q2);

//========================================================= Collins partcontent KT
void CollinsDistributionKt( struct PARTONCONTENT& fragmentation, PARAMETERS& Params, double z, double kt, double Q2);

//========================================================= Collins partcontent First moment
void CollinsDistributionFirstMoment( struct PARTONCONTENT& fragmentation,  PARAMETERS& Params, double z, double Q2);

//========================================================= Draw Collins partcontent
TCanvas* DrawCollinsDistribution( struct PARTONCONTENT& fragmentation,  PARAMETERS& Params, double Q2);



#endif //#ifndef __COLLINS_H_