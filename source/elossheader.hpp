#ifndef HEADERFILE_ELOSSHEADER
#define HEADERFILE_ELOSSHEADER

#include "linearinterpolation.hpp"
#include "grids.hpp"

//defining global variables to be used in all source files:

extern std::string collsys;				    		  	      //collision system
extern std::string sNN; 					    		  	  //collision energy
extern double nf;									          //effective number of flavours
extern std::string pName; 					    		  	  //particle name
extern std::string centrality;				    		  	  //centrality class
extern double xB;						    		  	      //xB value
extern double BCPP;									  	      //binary collision points percentage
extern size_t eventN;									  	  //number of events 
extern gridPoints Grids;				    		  	      //grid points
extern interpolationF dsdpti2; 				    		  	  //initial pT distribution
extern interpolationF LNorm, Ldndx, LColl;	    		  	  //interpolated L tables
extern double tau0;									  	      //thermalization time
extern double mgC, MC;					   			  	      //constant particle and gluon masses used for dA integrals
extern double TCollConst;				    		  	      //constant temperature used for Gauss filter integration
extern size_t phiGridN;				   					      //phi points number
extern std::vector<double> phiGridPts; 						  //defining vectors that store initial position points and angles
extern double TIMESTEP, TCRIT;							      //defining time step and critical temperature
extern double temp_tau0, temp_tauStep;
extern size_t temp_tauMax; 		 						      //temperature evolution tau grid parameters
extern double temp_x0,temp_xStep;
extern size_t temp_xMax; 		 							  //temperature evolution   x grid parameters
extern double temp_y0,temp_yStep;
extern size_t temp_yMax; 		 							  //temperature evolution   y grid parameters
extern std::vector<double> tempTauGrid, tempXGrid, tempYGrid; //temperature evolutions grids
extern int BCPSEED;										      //seed for generating initial position points

#endif