//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#include "unpolarised.h"
//#include <string>
#include <iostream>

using namespace std;

//========================================================= Unpolarized partcontent
void unpolarised(struct PARTONCONTENT& partcontent, double x, double Q2)
  // Returns the parton content at x,Q2
{
  std::vector<double> dens;


  for (int i=0; i<13; i++) dens.push_back(0.);

  int ISET = 1; //LO
  double UV, DV, US, DS, SS, GL;

  GRV98PA(ISET, x, Q2, UV, DV, US, DS, SS, GL);

  dens[8] = UV+US;
  dens[7] = DV+DS;
  dens[4] = US;
  dens[5] = DS;
  dens[9] = SS;
  dens[3] = SS;
  dens[10]= 0.;
  dens[2] = 0.;
  dens[11]= 0.;
  dens[1] = 0.;
  dens[12]= 0.;
  dens[0] = 0.;
  dens[6] = GL;
  
  //Assign parton densities PROTON:
  partcontent.up           =  dens[8] / x; // 
  partcontent.down         =  dens[7] / x; // 
  partcontent.anti_up      =  dens[4] / x; //
  partcontent.anti_down    =  dens[5] / x; //
  partcontent.strange      =  dens[9] / x; //
  partcontent.anti_strange =  dens[3] / x; //
  partcontent.charm        =  dens[10]/ x; //
  partcontent.anti_charm   =  dens[2] / x; //
  partcontent.bottom       =  dens[11]/ x; //
  partcontent.anti_bottom  =  dens[1] / x; //
  partcontent.top          =  dens[12]/ x; //
  partcontent.anti_top     =  dens[0] / x; //
  partcontent.glu          =  dens[6] / x; //
};


// definitions of DEUTRON
void Deutron( struct PARTONCONTENT& partcontent ){
    
    partcontent.up      += partcontent.down;
    partcontent.down     = partcontent.up;
    partcontent.anti_up   += partcontent.anti_down;
    partcontent.anti_down  = partcontent.anti_up;
    partcontent.strange     *= 2.;
    partcontent.anti_strange  *= 2.;
    partcontent.charm     *= 2.;
    partcontent.anti_charm  *= 2.;
    partcontent.glu     *= 2.;
    partcontent.bottom     *= 2.;
    partcontent.anti_bottom  *= 2.;
    partcontent.top     *= 2.;
    partcontent.anti_top  *= 2.;
    // IF WE USE VALENCE PARAMETRIZATION --> FORMULAS ARE THE SAME!

};

// definitions of NEUTRON
void Neutron( struct PARTONCONTENT& partcontent ){
      double neutron_pdf       = partcontent.up;
      partcontent.up    = partcontent.down;
      partcontent.down  = neutron_pdf;
      neutron_pdf       = partcontent.anti_up;
      partcontent.anti_up = partcontent.anti_down;
      partcontent.anti_down  = neutron_pdf;
    // IF WE USE VALENCE PARAMETRIZATION --> FORMULAS ARE THE SAME!
 };

// definitions of ANTIPROTON
void Antiproton( struct PARTONCONTENT& partcontent ){
      double antiproton_pdf       = partcontent.up;
      partcontent.up    = partcontent.anti_up;
      partcontent.anti_up = antiproton_pdf;

      antiproton_pdf    = partcontent.down;
      partcontent.down  = partcontent.anti_down;
      partcontent.anti_down= antiproton_pdf;

      antiproton_pdf    = partcontent.strange;
      partcontent.strange   = partcontent.anti_strange;
      partcontent.anti_strange  = antiproton_pdf;
}

double getTargetMass(std::string target)
{
  if (target == "proton")
  {
    return mpr;
  }
  else  if (target == "neutron")
  {
    return mneutron;
  }
  else  if (target == "deuteron")
  {
    return mdeuteron;
  }
  else  if (target == "antiproton")
  {
    return mpr;
  }
  else
  {
    std::cerr << "Target " << target << " is not supported." << std::endl;  
    std::cerr << "Use proton, neutron, deuteron, antiproton." << std::endl;   
    return mpr;
  }
  
}


enum TARGET_TYPE getTarget(std::string target)
{
  enum TARGET_TYPE target_type;

  if (target == "proton")
  {
    target_type = PROTON;
  }
  else  if (target == "neutron")
  {
    target_type = NEUTRON;
  }
  else  if (target == "deuteron")
  {
    target_type = DEUTERON;
  }
  else  if (target == "antiproton")
  {
    target_type = ANTIPROTON;
  }
  else
  {
    std::cerr << "Target " << target << " is not supported." << std::endl; 
    std::cerr << "Use proton, neutron, deuteron, antiproton." << std::endl;   
    target_type = PROTON;
  }

  return target_type;
  
}

//========================================================= Unpolarised partcontent
void UnpolarisedDistribution(std::string target, struct PARTONCONTENT& partcontent, double x, double Q2)
  // Returns the parton content at x,Q2
{

  unpolarised(partcontent, x, Q2);

  switch( getTarget(target) ){
  case PROTON:
    break;
  case ANTIPROTON:
    Antiproton(partcontent);
    break;
  case DEUTERON:
    Deutron(partcontent);
    break;
  case NEUTRON:
    Neutron(partcontent);
    break;
  default:
    break;
  }
  
  
};


//========================================================= fragmentation functions for pi+
void Fragmentation(struct PARTONCONTENT& fragmentation, double z, double Q2)
  // Returns fragmentation content at z,Q2
{
  double dff[11];

  dff[0] = 0.; 
  dff[1] = 0.; 
  dff[2] = 0.;
  dff[3] = 0.;
  dff[4] = 0.; 
  dff[5] = 0.;
  dff[6] = 0.; 
  dff[7] = 0.;
  dff[8] = 0.; 
  dff[9] = 0.; 
  dff[10]= 0.;



  int positive    = 1;
  //int negative    = 2;
  //int neutral     = 0;
  

   // ffset   1,2,3 means K, KKP, BFGW
   // fforder 0,1 is LO, NLO(MSbar)
   // ihadron 1,2,3,4,5 is pi,K,h,p,n
   // icp     1,2,3 chooses between particle, anti-particle or sum of both 
   // ipi     1,2,3 is a flag for BFGW; inactive for K, KKP  
   // icharge 0,1,2,3 is 0,+,-,+&-   
   // int ffset  = 1; // Kretzer
  int ffset  = 4; // DSS
  int fforder= 0; // LO 
  int ihadron= 1; // PION
  int icp    = 1; //fragmentation.icp;
  int ipi    = 1;  
  int icharge= positive;


  
  DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi);


  fragmentation.up            = dff[1+5];
  fragmentation.down          = dff[2+5];
  fragmentation.anti_up       = dff[-1+5];
  fragmentation.anti_down     = dff[-2+5];
  fragmentation.strange       = dff[3+5];
  fragmentation.anti_strange  = dff[-3+5];
  fragmentation.charm         = dff[4+5];
  fragmentation.anti_charm    = dff[-4+5];
  fragmentation.glu           = dff[0+5];
  fragmentation.bottom        = dff[5+5];
  fragmentation.anti_bottom   = dff[-5+5];
  fragmentation.top           = 0.;
  fragmentation.anti_top      = 0.;
       
}


enum HADRON_TYPE getHadron(std::string hadron)
{
  enum HADRON_TYPE hadron_type;
  if (hadron == "pi+")
  {
    return hadron_type = PIPLUS;
  }
  else  if (hadron == "pi-")
  {
    return hadron_type = PIMINUS;
  }
  else  if (hadron == "pi0")
  {
    return hadron_type =PIZERO;
  } 
  else if (hadron == "k+")
  {
    return hadron_type = KPLUS;
  }
  else  if (hadron == "k-")
  {
    return hadron_type =KMINUS;
  }
  else  if (hadron == "k0")
  {
    return hadron_type = KZERO;
  } 
  else if (hadron == "h+")
  {
    return hadron_type =HPLUS;
  }
  else  if (hadron == "h-")
  {
    return hadron_type = HMINUS;
  }
  else  if (hadron == "h0")
  {
    return hadron_type = HZERO;
  } 
  else
  {
    std::cerr << "Hadron " << hadron << " is not supported." << std::endl;
    std::cerr << "Use pi+,pi-,pi0, k+,k-, k0, h+,h0, h0" << std::endl;    
    return hadron_type = PIPLUS;
  }
  
}


//========================================================= fragmentation functions for all species
void Fragmentation(std::string hadron, struct PARTONCONTENT& fragmentation, double z, double Q2)
  // Returns fragmentation content at z,Q2
{
 double dff[11];

  dff[0] = 0.; 
  dff[1] = 0.; 
  dff[2] = 0.;
  dff[3] = 0.;
  dff[4] = 0.; 
  dff[5] = 0.;
  dff[6] = 0.; 
  dff[7] = 0.;
  dff[8] = 0.; 
  dff[9] = 0.; 
  dff[10]= 0.;



  int positive    = 1;
  int negative    = 2;
  int neutral     = 0;
  

   // ffset   1,2,3 means K, KKP, BFGW
   // fforder 0,1 is LO, NLO(MSbar)
   // ihadron 1,2,3,4,5 is pi,K,h,p,n
   // icp     1,2,3 chooses between particle, anti-particle or sum of both 
   // ipi     1,2,3 is a flag for BFGW; inactive for K, KKP  
   // icharge 0,1,2,3 is 0,+,-,+&-   
   // int ffset  = 1; // Kretzer
  int ffset  = 4; // DSS
  int fforder= 0; // LO 
  int ihadron= 1; // PION
  int icp    = 1; //fragmentation.icp;
  int ipi    = 1;  
  int icharge= positive;

  switch( getHadron(hadron) ){
  case PIPLUS:
  case PIMINUS:
  case PIZERO:
    ihadron = 1;
    break;
  case KPLUS:
  case KMINUS:
  case KZERO:
    ihadron = 2;
    break;
  case HPLUS:
  case HMINUS:
  case HZERO:
    ihadron = 3;
    break;
  default:
    break;        
  }

  switch( getHadron(hadron) ){
  case PIPLUS:
  case KPLUS:
  case HPLUS:
    icharge= positive;
    break;
  case PIMINUS:
  case KMINUS:
  case HMINUS:
    icharge= negative;
    break;
  case PIZERO:
  case HZERO:
  case KZERO:
    icharge= neutral;
    break;
  default:
    break;        
  }
  
  DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi);

   if( ihadron == 2 && icharge == 0 )
   { // Our assumptions
    // for K0_S = 1/sqrt{2} {d {bar s} + {bar d} s) !!!!
    
    double SB_Kp, DB_Kp, UB_Kp, GL_Kp, U_Kp, D_Kp, S_Kp,
      SB_Km, DB_Km, UB_Km, GL_Km, U_Km, D_Km, S_Km; 


    icharge = positive;
    DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi); // K^+ is produced!!!
  

    SB_Kp  = dff[5-3];
    DB_Kp  = dff[5-2];
    UB_Kp  = dff[5-1];
    GL_Kp  = dff[5+0];
    U_Kp   = dff[5+1];
    D_Kp   = dff[5+2];
    S_Kp   = dff[5+3];

    icharge = negative;
    DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi); // K^- is produced!!!

    SB_Km  = dff[5-3];
    DB_Km  = dff[5-2];
    UB_Km  = dff[5-1];
    GL_Km  = dff[5+0];
    U_Km   = dff[5+1];
    D_Km   = dff[5+2];
    S_Km   = dff[5+3];


    dff[5-5] = 0.  ;
    dff[5-4] = 0.  ;
    dff[5-3] = 0.5 * ( SB_Kp + SB_Km ); // D_{s-bar}^{K0} =   1/2[ D_{s-bar}^{K+} + D_{s-bar}^{K-}] ] 
    dff[5-2] = 0.5 * ( UB_Kp + UB_Km ); // D_{d-bar}^{K0} =  1/2[ D_{u-bar}^{K+} + D_{u-bar}^{K-}] ] ???
    dff[5-1] = 0.5 * ( DB_Kp + D_Km );  // D_{u-bar}^{K0} =  1/2[ D_{d-bar}^{K+} + D_{d}^{K-}] ]     ???
    dff[5+0] = 0.5 * ( GL_Kp + GL_Km );
    dff[5+1] = 0.5 * ( D_Kp  + DB_Km  );// D_{u}^{K0} =   1/2[ D_{d}^{K+} + D_{d-bar}^{K-}] ] 
    dff[5+2] = 0.5 * ( U_Kp  + U_Km  ); // D_{d}^{K0} =   1/2[ D_{u}^{K+} + D_{u-bar}^{K-}] ]  
    dff[5+3] = 0.5 * ( S_Kp  + S_Km  ); // D_{s}^{K0} =   1/2[ D_{s}^{K+} + D_{s}^{K-}] ]
    dff[5+4] = 0.  ;
    dff[5+5] = 0.  ;


    icharge = neutral;  // STATUS QUO IS RESTORED.
   }



  fragmentation.up            = dff[1+5];
  fragmentation.down          = dff[2+5];
  fragmentation.anti_up       = dff[-1+5];
  fragmentation.anti_down     = dff[-2+5];
  fragmentation.strange       = dff[3+5];
  fragmentation.anti_strange  = dff[-3+5];
  fragmentation.charm         = dff[4+5];
  fragmentation.anti_charm    = dff[-4+5];
  fragmentation.glu           = dff[0+5];
  fragmentation.bottom        = dff[5+5];
  fragmentation.anti_bottom   = dff[-5+5];
       

}