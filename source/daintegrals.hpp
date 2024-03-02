#ifndef HEADERFILE_DAHEADER
#define HEADERFILE_DAHEADER

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FdA integrals definitions:
void FdAHaltonSeqInit(size_t FdAMaxPts);										       //generates Halton sequences for FdA integrals
double dAp410(double ph, interpolationF &normint);								       //dAp410 integral definition
double FdA411(double ph, double dp, interpolationF &normint, interpolationF &dndxint); //FdA411 integral definition
double FdA412(double ph, double dp, interpolationF &normint, interpolationF &dndxint); //FdA412 integral definition
double FdA413(double ph, double dp, interpolationF &normint, interpolationF &dndxint); //FdA413 integral definition
double FdA414(double ph, double dp, interpolationF &normint, interpolationF &dndxint); //FdA414 integral definition
double FdA415(double ph, double dp, interpolationF &normint, interpolationF &dndxint); //FdA415 integral definition

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//dA integrals definitions:
void dAHaltonSeqInit(size_t dAMaxPts); 							           //generates Halton sequences for dA integrals
double dA410(double ph, interpolationF &normint); 					       //dA410 integral definition
double dA411(double ph, interpolationF &normint, interpolationF &dndxint); //dA411 integral definition
double dA412(double ph, interpolationF &normint, interpolationF &dndxint); //dA412 integral definition
double dA413(double ph, interpolationF &normint, interpolationF &dndxint); //dA413 integral definition
double dA414(double ph, interpolationF &normint, interpolationF &dndxint); //dA414 integral definition
double dA415(double ph, interpolationF &normint, interpolationF &dndxint); //dA415 integral definition
double dA416(double ph, interpolationF &normint, interpolationF &dndxint); //dA416 integral definition
double dA417(double ph, interpolationF &normint, interpolationF &dndxint); //dA417 integral definition

#endif