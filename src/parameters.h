//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#ifndef  __PARAMETERS_H_
#define  __PARAMETERS_H_ 

#include <vector>
#include <string>

// Parameters class
class VALUES {
 private:
 public:
  VALUES() {};
  ~VALUES() {};

  // Down quark
  double a_down;
  double b_down;
  double n_down; 
  double m2_down; // kt dependence


  // Up quark
  double a_up;
  double b_up;
  double n_up;
  double m2_up; // kt dependence

  // anti Down quark
  double a_anti_down;
  double b_anti_down;
  double n_anti_down; 
  double m2_anti_down; // kt dependence


  // anti Up quark
  double a_anti_up;
  double b_anti_up;
  double n_anti_up;
  double m2_anti_up; // kt dependence

  // Strange quark
  double a_strange;
  double b_strange;
  double n_strange; 
  double m2_strange; 


  // anti Strange quark
  double a_anti_strange;
  double b_anti_strange;
  double n_anti_strange;
  double m2_anti_strange;
};




const int set_number = 200;


class SETS{ // organize all sets here params, errors, sets of params
 private:
  VALUES main_parameters;
  std::string FILE_READ; // file to read parameters
  std::string SETS_FILE_READ; // file to read set parameters
  std::string FILE_WRITE; // file to write parameters
  std::string SETS_FILE_WRITE; // file to read set parameters

 public:
  SETS() {}; // constructor...
  ~SETS() {}; // deconstructor...

  VALUES parameters;
  VALUES errors;

  VALUES set[set_number];
  void SetReadFile(std::string file);
  void SetSetsReadFile(std::string file);
  void SetWriteFile(std::string file);
  void SetSetsWriteFile(std::string file);
  void ReadParameters();
  void ReadSets();
  void WriteParameters();
  void WriteSets();
  void UseMainSet(void);
  void UseSet(int k);

  std::vector<double> MinMax( double (*f)(double x), double x ); // evaluate min and max of a function f(x) on the parameters sets


};


// HERE WE DEFINE ALL PARAMETERS...
class PARAMETERS{
 private:
  // All params...
  double kt2_average;
  double ptq2_average;


 public:
  PARAMETERS(); // constructor...
  //Get params...
  double GetKt2Average( void ) const { return kt2_average;};
  double GetPtq2Average( void ) const { return ptq2_average;};


  SETS Sivers;
  SETS Transversity;
  SETS Collins;
};

#endif // #ifndef  __PARAMETERS_H_
