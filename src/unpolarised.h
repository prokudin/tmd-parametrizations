//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef __UNPOLARISED_H_
#define __UNPOLARISED_H_
#include "cfortran.h"
//#include <string>
#include <vector>
#include <iostream>
#include <cmath>

/* * */
/*       SUBROUTINE GRV98PA (ISET, X, Q2, UV, DV, US, DS, SS, GL) */
/* ********************************************************************* */
/* *                                                                   * */
/* *   THE PARTON ROUTINE.                                             * */
/* *                                     __                            * */
/* *   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         * */
/* *            X  =  Bjorken-x        (between  1.E-9 and 1.)         * */
/* *            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         * */
/* *                                                                   * */
/* *   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       * */
/* *            DS = d(bar),  SS = s = s(bar),  GL = gluon.            * */
/* *            Always x times the distribution is returned.           * */
/* *                                                                   * */
/* *   COMMON:  The main program or the calling routine has to have    * */
/* *            a common block  COMMON / INTINIP / IINIP , and the     * */
/* *            integer variable  IINIP  has always to be zero when    * */
/* *            GRV98PA is called for the first time or when  ISET     * */
/* *            has been changed.                                      * */
/* *                                                                   * */
/* *   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   * */
/* *            (1+1809 lines with 6 columns, 4 significant figures)   * */
/* *                                                                   * */
/* *******************************************************i************* */
/* * */
PROTOCCALLSFSUB9(GRV98PA,grv98pa,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define GRV98PA(ISET, x, Q2, UV, DV, US, DS, SS, GL) \
CCALLSFSUB9(GRV98PA,grv98pa,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, x, Q2, UV, DV, US, DS, SS, GL)


// FRAGMENTATION LIB
/* ************************************************************************ */
/* ************************************************************************ */
/* * Stefan Kretzer (kretzer@pa.msu.edu) Nov 2001:                        * */
/* *                                                                      * */
/* * Library of FF parametrizations updated to be posted on               * */
/* * http://www.pv.infn.it/~radici/FFdatabase/                            * */
/* * site maintained by Marco Radici (radici@pv.infn.it)                  * */
/* * and Rainer Jakob (rainer@theorie.physik.uni-wuppertal.de)            * */
/* * within the EU network:                                               * */
/* * "Hadronic Physics with High Energy Electromagnetic Probes"           * */
/* *                                                                      * */
/* * Also posted on that site are the individual FF sets included in this * */
/* * library as distributed by the corresponding authors.                 * */
/* ************************************************************************ */
/* * K, KKP and BFGW fragmentation functions for                          * */
/* * iparton -> ihadron                                                   * */
/* *                                                                      * */
/* * Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        * */
/* * "dff(iparton)" of dimension dff(-5:5)                                * */
/* *                                                                      * */
/* * iparton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar      * */
/* *                                                                      * */
/* * [Q2] = GeV^2                                                         * */
/* *                                                                      * */
/* * ffset = 1,2,3 means K, KKP, BFGW                                     * */
/* ************************************************************************ */
/* *        ALEXEI PROKUDIN prokudin@to.infn.it                           * */
/* *        making the same interface FF fDSS as for DLIB                 * */
/* ************************************************************************ */
/* * ffset = 4 means                                                      * */
/* *        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                     * */
/* *  D.de Florian, R.Sassot, M.Stratmann   hep-ph/0703242)               * */
/* *    Phys.Rev.D.75:114010,2007                                         * */
/* *                                                                      * */
/* ************************************************************************ */
/* *        ALEXEI PROKUDIN prokudin@to.infn.it                           * */
/* *        making the same interface FF AKK as for DLIB                  * */
/* *                    26/03/2008                                        * */
/* * ffset = 5 means                                                      * */
/* c---------------------------------------------------------------------- */
/* * */
/* c     AKK ROUTINES 2008 */
/* * */
/* c----------------------------------------------------------------------* */
/* * */
/* *      SUBROUTINE AKK(IH,Z,Q,DH)* */
/* * */
/* * AKK Update: Improvements from New Theoretical Input and Experimental Data. */
/* * S. Albino, B.A. Kniehl, G. Kramer . Mar 2008.  */
/* * e-Print: arXiv:0803.2768 [hep-ph] */
/* * Details are to be found in the corresponding references:             * */
/* * S.Kretzer, Phys.Rev.D62, 054001 (2000)                               * */
/* * B.A.Kniehl, G.Kramer, B.Potter, Nucl.Phys.B582, 514 (2000)           * */
/* * L.Bourhis, M.Fontannaz, J.P.Guillet, M.Werlen,Eur.Phys.J.C19,89(2001)* */
/* *                                                                      * */
/* * fforder=0,1 is LO, NLO(MSbar)                                        * */
/* * no LO set for BFGW                                                   * */
/* *                                                                      * */
/* * ihadron=1,2,3,4,5 is pi,K,h,p,n                                      * */
/* * no pi,K sets for BFGW                                                * */
/* * no p,n sets for K, BFGW                                              * */
/* * no n set for fDSS                                                    * AP */
/* *                                                                      * */
/* * icharge=0,1,2,3 is 0,+,-,+&-                                         * */
/* * note: 3 = +&- = charge sum (NOT average) in this library             * */
/* * for the neutral (icharge=0) particles K^0, n:                        * */
/* * icp = 1,2,3 chooses between particle, anti-particle or sum of both   * */
/* * icp is inactive for pi,h,p and charged Kaons                         * */
/* *                                                                      * */
/* * Following KKP, FFs into {pi^0;(anti-)K^0} and into {pi^+/-;K^+/-}    * */
/* * are related by isospin.                                              * */
/* *                                                                      * */
/* * ipi = 1,2,3 is a flag for BFGW; inactive for K, KKP:                 * */
/* * ipi = 1: best fit     (formula (8))  of BFGW                         * */
/* * ipi = 2: large Ng set (formula (9))  of BFGW                         * */
/* * ipi = 3: low Ng set   (formula (10)) of BFGW                         * */
/* ************************************************************************ */
/* * If options (ihadron, icharge etc.) are chosen which do not exist for * */
/* * a given parametrization (ffset) or if a flag is chosen outside its   * */
/* * range as defined above then the code returns a corresponding warning * */
/* * and stops.                                                           * */
/* ************************************************************************ */
/* ************************************************************************ */
PROTOCCALLSFSUB9(DLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT)
#define DLIB(Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI) \
  CCALLSFSUB9(MYDLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT,\
	      Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI)


// Partonic content DIS
struct PARTONCONTENT {
  double up,down,anti_up,anti_down,strange,anti_strange,charm,anti_charm,bottom,anti_bottom,top,anti_top,glu;
};


const  double eu               =  2./3.; // Up quark charge
const  double ed               = -1./3.; // Down quark charge
const  double es               = -1./3.; // Strange quark charge
const double PI                = 3.1415926543; // Pi
const double mpr               = 0.93827203 ; // the proton mass
const double mdeuteron          = 1.875613/2.   ; // the deutron mass/2. !!!!!!!!!!!
const double mneutron          = 0.93956536;    // the neutron mass
const double mpion             = 0.13956995;	// Mass of pion
const double mkaon             = 0.493667;	// Mass of kaon+-
const  double alpha_em         =  1./137.035999679; // aplha_em0
// Fragmentation
// ffset   1,2,3 means K, KKP, BFGW
// fforder 0,1 is LO, NLO(MSbar)
// ihadron 1,2,3,4,5 is pi,K,h,p,n
// icp     1,2,3 chooses between particle, anti-particle or sum of both 
// ipi     1,2,3 is a flag for BFGW; inactive for K, KKP  
// icharge 0,1,2,3 is 0,+,-,+&-
const int pion       = 1;
const int hadron     = 2;
const int kaon       = 3;
const int proton     = 4;
const int deutron    = 5;
const int neutron    = 6;
const int antiproton = 7;

const double positive    = +1.;
const double negative    = -1.;
const double neutral     =  0.;

/// Possible targets
enum TARGET_TYPE 
{
  PROTON  = 1,
  DEUTERON,
  NEUTRON,
  ANTIPROTON
};

/// Possible hadrons
enum HADRON_TYPE 
{
  PIPLUS = 1,
  PIMINUS,
  PIZERO,
  KPLUS,
  KMINUS,
  KZERO,
  HPLUS,
  HMINUS,
  HZERO
};

void Deutron( struct PARTONCONTENT& partcontent );
void Neutron( struct PARTONCONTENT& partcontent );
void Antiproton( struct PARTONCONTENT& partcontent );
enum TARGET_TYPE getTarget(std::string target);
double getTargetMass(std::string target);
enum HADRON_TYPE getHadron(std::string hadron);


//========================================================= Unpolarized partcontent
void unpolarised(struct PARTONCONTENT& partcontent, double x, double Q2);

//========================================================= Unpolarised partcontent for all targets
void UnpolarisedDistribution(std::string target, struct PARTONCONTENT& partcontent, double x, double Q2);

//========================================================= fragmentation functions for pi+
void Fragmentation(struct PARTONCONTENT& fragmentation, double z, double Q2);

//========================================================= fragmentation functions for all hadrons
void Fragmentation(std::string hadron, struct PARTONCONTENT& fragmentation, double z, double Q2);

#endif //#ifndef __UNPOLARISED_H_