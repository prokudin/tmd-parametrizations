#ifndef __PARTICLE_H
#define __PARTICLE_H
//class particle

// --- This is for root use in C++
#include "TVector3.h"
#include "constants.h"


// Switch to longitudinal, transverse polarization......
enum POLARIZATION {
  unp = 1,
  longitudinal,
  transverse,
  NA
};


class PARTICLE {
 private:
  char     Name[32];
  double   Charge;
  double   Spin;
  double   Mass;
  TVector3 Momentum;
  TVector3 PolarizationVector;
  enum POLARIZATION Polarization;

 public:
  PARTICLE(); 
  virtual ~PARTICLE() {}; 
  PARTICLE(char* Name, double Charge, double Spin, double Mass);

  PARTICLE& operator=(const PARTICLE& P) {
    strcpy(Name, P.Name);
    Charge = P.Charge;
    Spin = P.Spin;
    Mass = P.Mass;
    Momentum = P.Momentum;
    PolarizationVector = P.PolarizationVector;
    Polarization = P.Polarization;
  }



  virtual void SetCharge(double _Charge) { Charge = _Charge;}; 
  virtual void  SetSpin(double _Spin) { Spin = _Spin;}; 
  virtual void  SetMass(double _Mass) { Mass = _Mass;}; 
  virtual void  SetName(char* _Name) { strcpy(Name, _Name);}; 



  double GetCharge() const {return Charge;}; // these are to read protected values...
  double GetSpin() const {return Spin;}; 
  double GetMass() const {return Mass;}; 

  enum POLARIZATION GetPolarization() const {return Polarization;};
  void SetPolarization( enum POLARIZATION p ) { 
    if ( GetSpin() == 0 ) Polarization = NA; 
    else Polarization = p;
  };

  virtual char*  GetName() {return Name;}
  virtual TVector3 GetMomentum() const {return Momentum;};
  virtual TVector3 GetPolarizationVector() const {return PolarizationVector;};

  //  virtual double GetEnergy() const { return mod3(Mass, Momentum);}

  virtual double GetEnergy() const { return sqrt(Mass*Mass + Momentum*Momentum);}

  void SetMomentum( const TVector3&  p) { Momentum = p; };
  void SetPolarizationVector( const TVector3& s) { PolarizationVector = s; };

   virtual void print( std::ostream& kuda ); 
//  virtual void print();
}; 

#endif // #ifndef __PARTICLE_H
