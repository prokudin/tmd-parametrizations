#include "particle.h"
#include <iostream>

using namespace std;


// // Constructor
PARTICLE::PARTICLE(){
  strcpy( this->Name,  "default");
  this->Charge     = 0.;
  this->Spin       = 0.;
  this->Mass       = 0.;
  this->Momentum   = TVector3(0.,0.,0.);
  this->PolarizationVector   = TVector3(0.,0.,0.);
};

// // Constructor
PARTICLE::PARTICLE(char* _Name, double _Charge = 0, double _Spin = 0, double _Mass = 0) :
  Charge(_Charge), Spin(_Spin), Mass(_Mass) {
  strcpy(Name, _Name);
};


void PARTICLE::print( ostream& kuda){
    kuda << " *************************************" << endl;
    kuda << " *           PARTICLE                *" << endl;
    kuda << " *  Name:    " << GetName() << "                  *"<< endl;
    kuda << " *  Charge:  " << GetCharge() << "                       *"<< endl;
    kuda << " *  Spin:    " << GetSpin() << "                       *"<< endl;
    kuda << " *  Mass:    " << GetMass() << "                       *"<< endl;
    kuda << " *  Momentum: ";     Momentum.Print();
    kuda << " *  PolarizationVector: ";     PolarizationVector.Print();
    kuda << " *  Polarization: ";
    switch( GetPolarization() ){
    case transverse:
      kuda << " transverse " << "                       *"<< endl;
      break;
    case longitudinal:
      kuda << " longitudinal " << "                       *" <<endl;
      break;
    case unp:
      kuda << " unpolarised " << "                       *"<< endl;
      break;
    default:
      kuda << " N/A         " << "                       *"<< endl;
      break;
    }
    kuda << " *************************************" << endl;
};
