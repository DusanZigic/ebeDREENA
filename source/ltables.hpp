#ifndef HEADERFILE_LTABLESHEADER
#define HEADERFILE_LTABLESHEADER

#include "grids.hpp"

//defining global variables to be used in all source files:

extern std::string LT_sNN;    //particle name
extern double LT_nf;          //effective number of flavours
extern std::string LT_pName;  //particle name
extern double LT_xB;          //xB value
extern gridPoints LT_Grids;   //grids
extern size_t LdndxMaxPoints; //maximal number of points for Ldndx integration
extern size_t LCollMaxPoints; //maximal number of points for collisional integration

#endif