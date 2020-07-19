Author: Alexei Prokudin, 
prokudin@jlab.org

Let me know if any bug is present / improvement of the code is welcome.


***************
Installation:

```make all```


***************
Cleaning:
```make clean```


***************
Usage:

```test.exe``` 

produces three plots on 

Collins FF   arXiv:0812.4366
Transversity arXiv:0812.4366
Sivers       Eur.Phys.J.A39:89-100,2009


```stfunctions.exe```

prints values for structure functions

***************
Please cite **Anselmino:2008jk**
https://inspirehep.net/literature/806038

and **Anselmino:2008sga** 
https://inspirehep.net/literature/786122

***************
Usage of the structure functions:

*************** 
Sivers SF
***********************

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

*************** 
Collins SF
***********************
 
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

*************** 
Unpolarised SF
***********************

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
 
