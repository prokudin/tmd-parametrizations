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

//**************************************************************************************
//**************************          SIVERS      **************************************
//========================================================= Sivers effect...x dependence
double sivers_x_dependence(double x, double a, double b, double n);

//========================================================= Sivers effect...
void sivers(struct PARTONCONTENT& partcontent, PARAMETERS Params, double x );

//========================================================= Unpolarized partcontent
void unpolarised(struct PARTONCONTENT& partcontent, double x, double Q2);

//========================================================= Sivers partcontent
void SiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= Unpolarised partcontent
void UnpolarisedDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= Sivers partcontent KT
void SiversDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double kt, double Q2);

//========================================================= Sivers partcontent First moment
void SiversDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= Draw Sivers partcontent
TCanvas* DrawSiversDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double Q2);







//**************************************************************************************
//**************************          COLLINS     **************************************
//========================================================= Collins effect...z dependence
double collins_z_dependence(double z, double a, double b, double n);

//========================================================= Fragmentation PI^+
void Fragmentation(struct PARTONCONTENT& fragmentation, double z, double Q2);

//========================================================= Collins effect PI^+...
void collins(struct PARTONCONTENT& fragmentation, PARAMETERS Params, double z);

//========================================================= Collins partcontent
void CollinsDistribution( struct PARTONCONTENT& fragmentation, PARAMETERS Params, double z, double Q2);

//========================================================= Collins partcontent KT
void CollinsDistributionKt( struct PARTONCONTENT& fragmentation, PARAMETERS Params, double z, double kt, double Q2);

//========================================================= Collins partcontent First moment
void CollinsDistributionFirstMoment( struct PARTONCONTENT& fragmentation,  PARAMETERS Params, double z, double Q2);

//========================================================= Draw Collins partcontent
TCanvas* DrawCollinsDistribution( struct PARTONCONTENT& fragmentation,  PARAMETERS Params, double Q2);









//**************************************************************************************
//**************************          TRANSVERSITY      ********************************
//========================================================= Sivers effect...x dependence
double transversity_x_dependence(double x, double a, double b, double n);

//========================================================= Soffer Bound...
void SofferBound(struct PARTONCONTENT& partcontent, double x, double Q2 );

//========================================================= transversity
void transversity(struct PARTONCONTENT& partcontent, PARAMETERS Params, double x );

//========================================================= partcontent
void TransversityDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);


//========================================================= partcontent KT
void TransversityDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double kt, double Q2);

//========================================================= partcontent First moment
void TransversityDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

//========================================================= Draw partcontent
TCanvas* DrawTransversityDistribution( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double Q2);
