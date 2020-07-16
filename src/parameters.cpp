//
// Author: Alexei Prokudin <prokudin@jlab.org>
//
#include "parameters.h"

#include <iostream>
#include <fstream>

using namespace std;

PARAMETERS::PARAMETERS() {
  // init all parameters
  kt2_average  = 0.25;
  ptq2_average = 0.2;

  Sivers.SetReadFile("sivers_parameters.dat");
  Sivers.ReadParameters();

  Transversity.SetReadFile("transversity_parameters.dat");
  Transversity.ReadParameters();


  Collins.SetReadFile("collins_parameters.dat");
  Collins.ReadParameters();

};



//========================================================= parameters
void SETS::ReadParameters(){

   ifstream in;
   in.open(FILE_READ.c_str()); // fit 

   if(!in) {
     cerr << "PARAMETERS: Unable to open file with parameters! " << FILE_READ << endl;
   } 

   in >>  parameters.a_up >>  errors.a_up;
   in >>  parameters.b_up >>  errors.b_up;
   in >>  parameters.n_up >>  errors.n_up;
   in >>  parameters.m2_up>>  errors.m2_up;

   in >>  parameters.a_down     >>  errors.a_down;
   in >>  parameters.b_down     >>  errors.b_down;
   in >>  parameters.n_down     >>  errors.n_down;
   in >>  parameters.m2_down    >>  errors.m2_down;

   in >>  parameters.a_anti_up >>  errors.a_anti_up;
   in >>  parameters.b_anti_up >>  errors.b_anti_up;
   in >>  parameters.n_anti_up >>  errors.n_anti_up;
   in >>  parameters.m2_anti_up>>  errors.m2_anti_up;

   in >>  parameters.a_anti_down     >>  errors.a_anti_down;
   in >>  parameters.b_anti_down     >>  errors.b_anti_down;
   in >>  parameters.n_anti_down     >>  errors.n_anti_down;
   in >>  parameters.m2_anti_down    >>  errors.m2_anti_down;


   in >>  parameters.a_strange >>  errors.a_strange;
   in >>  parameters.b_strange >>  errors.b_strange;
   in >>  parameters.n_strange >>  errors.n_strange;
   in >>  parameters.m2_strange>>  errors.m2_strange;

   in >>  parameters.a_anti_strange >>  errors.a_anti_strange;
   in >>  parameters.b_anti_strange >>  errors.b_anti_strange;
   in >>  parameters.n_anti_strange >>  errors.n_anti_strange;
   in >>  parameters.m2_anti_strange>>  errors.m2_anti_strange;


   in.close();

   main_parameters = parameters; // remember the main set of parameters
};



//========================================================= parameters
void SETS::UseSet(int k){
   parameters =  set[k];
}

//========================================================= parameters
void SETS::UseMainSet(void){
   parameters = main_parameters;
}


//========================================================= Constructor
void SETS::SetReadFile(std::string file){
  FILE_READ = file;
};

//========================================================= Constructor
void SETS::SetSetsReadFile(std::string file){
  SETS_FILE_READ = file;
};


//========================================================= Constructor
void SETS::SetWriteFile(std::string file){
  FILE_WRITE = file;
};

//========================================================= Constructor
void SETS::SetSetsWriteFile(std::string file){
  SETS_FILE_WRITE = file;
};

//========================================================= parameters
void SETS::ReadSets(){

   ifstream in;
   in.open(SETS_FILE_READ.c_str()); // fit 

   if(!in) {
     cerr << "PARAMETERS: SETS: Unable to open file with parameters! " << SETS_FILE_READ << endl;
   } 

   for(int i = 0; i < set_number; i++){
   in >>  set[i].a_up;
   in >>  set[i].b_up;
   in >>  set[i].n_up;
   in >>  set[i].m2_up;

   in >>  set[i].a_down;
   in >>  set[i].b_down;
   in >>  set[i].n_down;
   in >>  set[i].m2_down;

   in >>  set[i].a_anti_up;
   in >>  set[i].b_anti_up;
   in >>  set[i].n_anti_up;
   in >>  set[i].m2_anti_up;

   in >>  set[i].a_anti_down;
   in >>  set[i].b_anti_down;
   in >>  set[i].n_anti_down;
   in >>  set[i].m2_anti_down;


   in >>  set[i].a_strange;
   in >>  set[i].b_strange;
   in >>  set[i].n_strange;
   in >>  set[i].m2_strange;

   in >>  set[i].a_anti_strange;
   in >>  set[i].b_anti_strange;
   in >>  set[i].n_anti_strange;
   in >>  set[i].m2_anti_strange;
   }

   in.close();
};



//========================================================= parameters
void SETS::WriteParameters(){

   ofstream out;
   out.open(FILE_WRITE.c_str()); // fit 

   if(!out) {
     cerr << "PARAMETERS: Unable to open file to write parameters! " << FILE_WRITE << endl;
   } else {
     cout << "PARAMETERS: Open " << FILE_WRITE << endl;
   }

   out <<  parameters.a_up << " " << errors.a_up << endl;
   out <<  parameters.b_up << " " <<  errors.b_up << endl;
   out <<  parameters.n_up << " " <<  errors.n_up << endl;
   out <<  parameters.m2_up<< " " <<  errors.m2_up<< endl;

   out <<  parameters.a_down     << " " <<  errors.a_down << endl;
   out <<  parameters.b_down     << " " <<  errors.b_down << endl;
   out <<  parameters.n_down     << " " <<  errors.n_down << endl;
   out <<  parameters.m2_down    << " " <<  errors.m2_down << endl;

   out <<  parameters.a_anti_up << " " <<  errors.a_anti_up << endl;
   out <<  parameters.b_anti_up << " " <<  errors.b_anti_up << endl;
   out <<  parameters.n_anti_up << " " <<  errors.n_anti_up << endl;
   out <<  parameters.m2_anti_up<< " " <<  errors.m2_anti_up << endl;

   out <<  parameters.a_anti_down     << " " <<  errors.a_anti_down << endl;
   out <<  parameters.b_anti_down     << " " <<  errors.b_anti_down << endl;
   out <<  parameters.n_anti_down     << " " <<  errors.n_anti_down << endl;
   out <<  parameters.m2_anti_down    << " " <<  errors.m2_anti_down << endl;


   out <<  parameters.a_strange << " " <<  errors.a_strange << endl;
   out <<  parameters.b_strange << " " <<  errors.b_strange << endl;
   out <<  parameters.n_strange << " " <<  errors.n_strange << endl;
   out <<  parameters.m2_strange<< " " <<  errors.m2_strange << endl;

   out <<  parameters.a_anti_strange << " " <<  errors.a_anti_strange << endl;
   out <<  parameters.b_anti_strange << " " <<  errors.b_anti_strange << endl;
   out <<  parameters.n_anti_strange << " " <<  errors.n_anti_strange << endl;
   out <<  parameters.m2_anti_strange<< " " <<  errors.m2_anti_strange << endl;


   out.close();

};

//========================================================= parameters
void SETS::WriteSets(){

   ofstream out;
   out.open(SETS_FILE_WRITE.c_str()); // fit 

   if(!out) {
     cerr << "PARAMETERS: SETS: Unable to open file to write parameters! " << SETS_FILE_WRITE << endl;
   } else {
     cout << "PARAMETERS: SETS: Open " << SETS_FILE_WRITE << endl;
   }

   for(int i = 0; i < set_number; i++){
     out <<  set[i].a_up << " " ;
     out <<  set[i].b_up << " " ;
     out <<  set[i].n_up << " " ;
     out <<  set[i].m2_up << " " ;
     
     out <<  set[i].a_down << " " ;
     out <<  set[i].b_down << " " ;
     out <<  set[i].n_down << " " ;
     out <<  set[i].m2_down << " " ;
     
     out <<  set[i].a_anti_up << " " ;
     out <<  set[i].b_anti_up << " " ;
     out <<  set[i].n_anti_up << " " ;
     out <<  set[i].m2_anti_up << " " ;
   
     out <<  set[i].a_anti_down << " " ;
     out <<  set[i].b_anti_down << " " ;
     out <<  set[i].n_anti_down << " " ;
     out <<  set[i].m2_anti_down << " " ;


     out <<  set[i].a_strange << " " ;
     out <<  set[i].b_strange << " " ;
     out <<  set[i].n_strange << " " ;
     out <<  set[i].m2_strange << " " ;

     out <<  set[i].a_anti_strange << " " ;
     out <<  set[i].b_anti_strange << " " ;
     out <<  set[i].n_anti_strange << " " ;
     out <<  set[i].m2_anti_strange;
     out << endl;
   }

   out.close();
};


// Evaluate min and max of a function f(x) on the parameters sets
std::vector<double> SETS::MinMax( double (*f)(double ), double x ){
  std::vector<double> r;
  double min, max;

   for(int k = 0; k < set_number; k++){

    SETS::UseSet(k); // we do loop over all sets...

    double sechenie = (*f)( x );

    cout << k << endl;

    if( k == 0 ){
      min = sechenie;
      max = sechenie;
    } else {
      if( (sechenie = (*f)( x ) ) > max ){
	      max = sechenie;
      }
      if(  sechenie < min ){
	      min = sechenie;
      }
    }
  }

  SETS::UseMainSet(); // main set is restored after the loop ...
  
  r.push_back(min);
  r.push_back(max);
  
  return r;
      
};
