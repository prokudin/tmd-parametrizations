//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#include "stfunctions.h"
#include "sivers.h"
#include "transversity.h"
#include "collins.h"

double TMD::Sum(PARTONCONTENT& Target, PARTONCONTENT& Produced)
{    
    return 
      eu*eu * Target.up     *   Produced.up +
      ed*ed * Target.down   *   Produced.down +
      eu*eu * Target.anti_up  *   Produced.anti_up +
      ed*ed * Target.anti_down*   Produced.anti_down +
      es*es * Target.strange    *   Produced.strange +
      es*es * Target.anti_strange *   Produced.anti_strange; 
} // sum_q e_q * e_q * q(x) * D(z)

TMD::TMD()
{
    PARAMETERS P;
    this->Params = P;
}

/// FUT^sin(Phi_h + Phi_S) structure function
/// Parameters: 
/// target = proton, neutron, deuteron, antiproton
/// hadron = pi+,pi-,pi0, k+,k-, k0, h+,h0, h0
/// S energy in GeV2
/// x
/// z
/// Q2 in GeV2
/// PhT in GeV
double TMD::FUTCollins(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT)
{
  double kt2_average  = Params.GetKt2Average();
  double ptq2_average = Params.GetPtq2Average();
  double Q = sqrt(Q2);

  PARTONCONTENT Target, Produced;

  double targetmass = getTargetMass(target);

  double y = Q2/((S-pow2(targetmass))*x);

  double M2collins;

  if( getHadron(hadron) == PIPLUS || getHadron(hadron) == PIMINUS || getHadron(hadron) == PIZERO ) 
  {
    M2collins = Params.Collins.parameters.m2_up;
  }
  else
  {
    std::cerr << "Hadron " << hadron << " is not available for Collins FF." << std::endl;
    return 0;
  }
  
  double collins_width = M2collins * ptq2_average/(ptq2_average + M2collins);
  double dptav_numerator =  collins_width +  pow2(z)*Params.Transversity.parameters.m2_up ;

// It is the cross section here
//    double coeff_coll = (1.-y)*collins_width*PhT/(PI*pow2(dptav_numerator))*
//      	exp( -PhT*PhT/dptav_numerator)/(Q*Q*Q*Q)*S*x* PhT * 2. * pow2(2. * PI)  * pow2( alpha_em ); 

// Structure Function:
  double coeff_coll = collins_width*PhT/(PI*pow2(dptav_numerator))*exp( -PhT*PhT/dptav_numerator); 

  coeff_coll *=  sqrt( 2.*exp(1.)/ M2collins ) *collins_width/ptq2_average;
  coeff_coll *=  0.5; // in Eq (20) there's sum over N(z) D(z) I use 2 N(z) D(z) in CollinsDistribution 

  TransversityDistribution( target, Target, Params, x, Q2); // Partcontent for distribution
  
  CollinsFragmentation( hadron, Produced,  Params, z, Q2); // Partcontent for fragmentation TODO not working properly yet...

  return  coeff_coll * Sum(Target, Produced); //  sum_q e_q * e_barq * q(x) * barq(x)
}

/// FUT^sin(Phi_h - Phi_S) structure function
/// Parameters: 
/// target = proton, neutron, deuteron, antiproton
/// hadron = pi+,pi-,pi0, k+,k-, k0, h+,h0, h0
/// S energy in GeV2
/// x
/// z
/// Q2 in GeV2
/// PhT in GeV
double TMD::FUTSivers(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT)
{
  double kt2_average  = Params.GetKt2Average();
  double ptq2_average = Params.GetPtq2Average();
  double Q = sqrt(Q2);

  PARTONCONTENT Target, Produced;

  double targetmass = getTargetMass(target);

  double y = Q2/((S-pow2(targetmass))*x);
  double dn1 = sqrt( 2.*exp(1.) )/sqrt(  Params.Sivers.parameters.m2_up );
  double dkt1 = kt2_average * Params.Sivers.parameters.m2_up / ( kt2_average + Params.Sivers.parameters.m2_up );
  double dpt1 = ptq2_average + pow2(z)*dkt1;

// It is the cross section here
//    double coeff_siv = 0.5 * pow2(dkt1) / kt2_average * z * PhT / dpt1 * dn1 *
//   	(1. + pow2(1.-y) ) * exp(-PhT*PhT/dpt1)/dpt1/(Q*Q*Q*Q)*S*x*PhT*  (1./(2.*PI)) * 2.* pow2(2. * PI)  * pow2( alpha_em ); 
// Structure Function:
  double coeff_siv = 0.5 * pow2(dkt1) / kt2_average * z * PhT / dpt1 * dn1 * exp(-PhT*PhT/dpt1)/dpt1; 

  SiversDistribution( target, Target,  Params, x, Q2); // Partcontent for distribution
  
  Fragmentation( hadron, Produced, z, Q2); // Partcontent for fragmentation

  return  coeff_siv * Sum(Target, Produced); //  sum_q e_q * e_barq * q(x) * barq(x)
}


/// FUU unpolarised structure function
/// Parameters: 
/// target = proton, neutron, deuteron, antiproton
/// hadron = pi+,pi-,pi0, k+,k-, k0, h+,h0, h0
/// S energy in GeV2
/// x
/// z
/// Q2 in GeV2
/// PhT in GeV
double TMD::FUU(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT)
{
  double kt2_average  = Params.GetKt2Average();
  double ptq2_average = Params.GetPtq2Average();
  double Q = sqrt(Q2);

  PARTONCONTENT Target, Produced;

  double targetmass = getTargetMass(target);

  double y = Q2/((S-pow2(targetmass))*x);
  double dptav= ptq2_average + pow2(z) * kt2_average;
// It is the cross section here
//  double coeff_unp = (1.+pow2(1.-y))*exp(-PhT*PhT/dptav)/dptav/(Q*Q*Q*Q)*S*x*PhT* (1./(2.*PI)) *2.* pow2(2. * PI)  * pow2( alpha_em ); // It is the cross section here
// Structure Function:
  double coeff_unp = exp(-PhT*PhT/dptav)/dptav; // It is the cross section here

  UnpolarisedDistribution( target, Target, x, Q2); // Partcontent for distribution
  
  Fragmentation( hadron, Produced, z, Q2); // Partcontent for fragmentation

  return  coeff_unp * Sum(Target, Produced); //  sum_q e_q * e_barq * q(x) * barq(x)
}