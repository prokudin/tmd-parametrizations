Author: Alexei Prokudin, 
prokudin@jlab.org

Let me know if any bug is present / improvement of the code is welcome.


***************
Installation:

unzip torino.zip
make all





***************
Usage:

test.exe 

produces three plots on 

Collins FF   arXiv:0812.4366
Transversity arXiv:0812.4366
Sivers       Eur.Phys.J.A39:89-100,2009





***************
Usage of the distributions:

*************** 
Sivers distridution
***********************

kt dependent: 

void SiversDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double kt, double Q2);

output is "partcontent" including all flavours 

struct PARTONCONTENT {
  double up,down,anti_up,anti_down,strange,anti_strange,charm,anti_charm,bottom,anti_bottom,top,anti_top,glu;
};

kt integrated (int d^2 kt) as described in various papers

void SiversDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

*************** 
Collins FF 
***********************
ATTENTION: produces quark --> pi^+ fragmentation

kt dependent 
void CollinsDistributionKt( struct PARTONCONTENT& fragmentation, PARAMETERS Params, double z, double kt, double Q2);

kt integrated (int d^2 kt)
void CollinsDistributionFirstMoment( struct PARTONCONTENT& fragmentation,  PARAMETERS Params, double z, double Q2);

*************** 
Transversity 
***********************
kt dependent
void TransversityDistributionKt( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double kt, double Q2);

kt integrated (int d^2 kt)
void TransversityDistributionFirstMoment( struct PARTONCONTENT& partcontent,  PARAMETERS Params, double x, double Q2);

