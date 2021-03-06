//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef __STFUNCTIONS_H_
#define __STFUNCTIONS_H_
#include <string>
#include "parameters.h"
#include "unpolarised.h"
#include "sivers.h"


inline double pow2(double a) {
  return a*a;
};

class TMD 
{
 private: 
    PARAMETERS Params;

 public:
    TMD();
    ~TMD(){};

    double Sum(PARTONCONTENT& Target, PARTONCONTENT& Produced);

    PARTONCONTENT SIDIS(PARTONCONTENT& Target, PARTONCONTENT& Produced);

    double SIDIS_SUM(PARTONCONTENT& out);

    PARTONCONTENT FUUparton(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);
    double FUU(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);

    PARTONCONTENT FUTSiversparton(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);
    double FUTSivers(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);

    PARTONCONTENT FUTCollinsparton(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);
    double FUTCollins(std::string & target, std::string & hadron, double S, double x, double z, double Q2, double PhT);

};

#endif // #ifndef __STFUNCTIONS_H_