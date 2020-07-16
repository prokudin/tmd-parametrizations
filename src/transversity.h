//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef __TRANSVERSITY_H_
#define __TRANSVERSITY_H_
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


/* ********************************************************************* */
/* *                                                                   * */
/* *                    TRANSVERSITY DENSITIES                         * */
/* *                                                                   * */
/* *   INPUT:   ISET = number of the parton set :                      * */
/* *              ISET = 1  NEXT-TO-LEADING ORDER  (MS-bar)            *  */
/* *                        (DATA FILE 'transmaxnlo.grid' UNIT=11, TO  * */
/* *                         BE DEFINED BY THE USER )                  *  */
/* *              ISET = 2  LEADING ORDER                              *  */
/* *                        (DATA FILE 'transmaxlo_new.grid' UNIT=22   * */
/* *                                                                   * */
/* *            X  = Bjorken-x       (between  1.E-4  and  1)          * */
/* *            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        * */
/* *             (for values outside the allowed range the program     * */
/* *              writes a warning and extrapolates to the x and       * */
/* *              Q2 values requested)                                 * */
/* *                                                                   * */
/* *   OUTPUT:  U = x * DELTA u                                        *  */
/* *            D = x * DELTA d                                        * */
/* *            UB = x * UBAR                                          * */
/* *            DB = x * DBAR                                          *    */
/* *            ST = x * DELTA STRANGE = x * DELTA STRANGE(BAR)        *      */
/* *                                                                   * */
/* *          (  For the parton distributions always x times           * */
/* *                   the distribution is returned   )                * */
/* *                                                                   * */
/* *                                                                   * */
/* *   COMMON:  The main program or the calling routine has to have    * */
/* *            a common block  COMMON / INTINI / IINI , and  IINI     * */
/* *            has always to be zero when PARPOL is called for the    * */
/* *            first time or when 'ISET' has been changed.            * */
/* *                                                                   * */
/* ********************************************************************* */
PROTOCCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  PARPOLT(ISET, X, Q2, U, D, UB, DB, ST) \
CCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, X, Q2, U, D, UB, DB, ST)




//**************************************************************************************
//**************************          TRANSVERSITY      ********************************
//========================================================= x dependence
double transversity_x_dependence(double x, double a, double b, double n);

//========================================================= Soffer Bound...
void SofferBound(struct PARTONCONTENT& partcontent, double x, double Q2 );

//========================================================= transversity
void transversity(struct PARTONCONTENT& partcontent, PARAMETERS Params, double x );

void TransversityDistribution(std::string target, struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= partcontent
void TransversityDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);


//========================================================= partcontent KT
void TransversityDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double kt, double Q2);

//========================================================= partcontent First moment
void TransversityDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= Draw partcontent
TCanvas* DrawTransversityDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double Q2);


#endif //#ifndef __TRANSVERSITY_H_