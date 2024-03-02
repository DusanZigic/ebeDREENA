#ifndef HEADERFILE_INPUTHEADER
#define HEADERFILE_INPUTHEADER

#include <string>

int GetInputs(int argc, char const *argv[], std::string &csys, std::string &snn, std::string &pname, std::string &cent, double &xb, size_t &eventn, double &bcpp, size_t &phiptsn, double &tstep, double &tcrit, int &bcpseed); //function that gets inputs for energy loss calculations
int GetInputs(int argc, char const *argv[],                    std::string &snn, std::string &pname,                    double &xb, size_t &LdndxMaxPts, size_t &LCollMaxPts, double &tcrit); 				    				//function that gets inputs for LTables calculations

#endif