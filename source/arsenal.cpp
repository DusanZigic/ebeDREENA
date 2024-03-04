#include "arsenal.hpp"
#include "elossheader.hpp"
#include "importexport.hpp"
#include "linearinterpolation.hpp"
#include "polyintegration.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <cfloat>

double productLog(double x)
{
	if (x == 0.0) {
		return 0.0;
	}

	double w0, w1;
	if (x > 0.0) {
		w0 = std::log(1.2 * x / std::log(2.4 * x / std::log1p(2.4 * x)));
	}
	else {
		double v = 1.4142135623730950488 * std::sqrt(1.0 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (true) {
		double e = std::exp(w0);
		double f = w0 * e - x;
		w1 = w0 - f / ((e * (w0 + 1.0) - (w0 + 2.0) * f / (w0 + w0 + 2.0)));
		if (std::abs(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}
	return w1;
}

double unitStep(double x) {
    return (x < 0.0) ? 0.0 : 1.0;
}
long double unitStep(long double x) {
    return (x < 0.0L) ? 0.0L : 1.0L;
}

double haltonSequence(int index, int base)
{
	double f = 1.0;
	double res = 0.0;

	while (index > 0) {
		f = f / static_cast<double>(base);
		res += f * static_cast<double>(index % base);
		index = index / base; // integer division
	}

	return res;
}

int generateTempGrid()
{
	if (loadTempGridParams() != 1) return -1;

	for (size_t iTau=0; iTau<temp_tauMax; iTau++) {
		for (size_t iX=0; iX<temp_xMax; iX++) {
			for (size_t iY=0; iY<temp_yMax; iY++) {
				tempTauGrid.push_back(temp_tau0 + iTau*temp_tauStep);
				  tempXGrid.push_back(  temp_x0 +   iX*temp_xStep);
				  tempYGrid.push_back(  temp_y0 +   iY*temp_yStep);
			}
		}		
	}

	return 1;
}

void generatePhiGridPts() {
	for (size_t iPhi=0; iPhi<phiGridN; iPhi++) phiGridPts.push_back(2.0*M_PI*static_cast<double>(iPhi)/static_cast<double>(phiGridN-1));
}

int BCPSEED; //seed for generating initial position points

int generateInitPosPoints(size_t event_id, std::vector<double> &xPoints, std::vector<double> &yPoints)
{
	std::vector<std::vector<double>> bcpts; if (loadBinCollPoints(event_id, bcpts) != 1) return -1;

	size_t bsptsNum = static_cast<size_t>(BCPP*bcpts.size());

	if (bsptsNum < 1) bsptsNum = 1;

	if (BCPSEED == 0) {
		std::random_device rd; auto rng = std::default_random_engine{rd()};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}
	else {
		auto rng = std::default_random_engine{static_cast<long unsigned int>(BCPSEED)};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}

	for (size_t iBCP=0; iBCP<bsptsNum; iBCP++) {
		xPoints.push_back(bcpts[iBCP][0]);
        yPoints.push_back(bcpts[iBCP][1]);
	}

	return 1;
}

static void generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab)
//function that generates sampling points for Gaussian integration
//qGTab, fGTab - vectors that store sampling point <- output
{	
	double sigmaNum = 3.5; //setting sigma
	double sigmaStep = 0.25; //setting step
	size_t GTabLen = 2 * static_cast<size_t>(sigmaNum / sigmaStep) + 1; //setting length of sampling points
	
	double GaussTabSum = 0.0; //setting normalization sum to zero
	
	for (size_t iG=0; iG<GTabLen; iG++) //calculating sampling points
	{
		qGTab.push_back(-1.0*sigmaNum + static_cast<double>(iG)*sigmaStep); //setting qGaussTab values
		fGTab.push_back(std::exp(-qGTab.back()*qGTab.back()/2.0));          //setting fGaussTab values
		GaussTabSum += fGTab.back();                                        //adding to normalization sum
	}
	
	for (size_t iG=0; iG<GTabLen; iG++)  //normalizing
	{
		fGTab[iG] /= GaussTabSum; //dividing fGaussTab values with total sum
	}
}

void gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2)
//function that performs Gauss filter integration - modefied pT integration algorithm
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF RadRelInt(Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (dsdpti2.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / dsdpti2.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF RadRelInt(Grids.RadPts(), Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();            //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (dsdpti2.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2.interpolation(pT) * GFSum);
			}
		}
	}
}

void gaussFilterIntegrate(const interpolationF &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2)
//function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
//dsdpti2lquark - light quark initial pT distribution      						  <- input
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF RadRelInt(Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (dsdpti2lquark.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF RadRelInt(Grids.RadPts(), Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (dsdpti2lquark.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
			}
		}
	}
}

void gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA)
//function that performs Gauss filter integration - default algorithm
//radiativeRAA  - raditive RAA 							   <- input
//collisionalEL - collisional energy loss				   <- input
//singRAA 		- RAA array after Gauss filter integration <- output
{
    interpolationF RadRelInt(Grids.RadPts(),   radiativeRAA);  //creating radiative RAA interpolated function
    interpolationF muCollInt(Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	double GFSum; //defining sum variable for Gauss filter

	double dpT; //defining pT and dpT variables

	double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value

	double sigmaColl; //defining variable for collisional sigma
	
	//Gauss filter
	for (const auto &pT : Grids.finPts())
	{
		GFSum = 0.0L;

		muCollCurrVal = muCollInt.interpolation(pT);

		sigmaColl = std::sqrt(2.0*TCollConst*muCollCurrVal);

		qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

		if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
			double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}		
		
		if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
			double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}
		
		//calculating Gauss filter
		for (size_t iG=0; iG<qGaussTab.size(); iG++)
		{
			dpT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
			GFSum += (dsdpti2.interpolation(pT + dpT)*RadRelInt.interpolation(pT + dpT)*(pT + dpT) / pT * fGaussTab[iG]);
		}

		singRAA.push_back(1.0 / dsdpti2.interpolation(pT) * GFSum);
	}
}

void calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp)
{
	interpolationF pathLenghDistInt(phiGridPts, pathLenghDist);
	avgPathLength.push_back(cubicIntegrate(phiGridPts, pathLenghDist)/2.0/M_PI);
	avgPathLength.push_back((pathLenghDistInt.interpolation(phiGridPts.front()) + pathLenghDistInt.interpolation(phiGridPts.back()))/2.0);
	avgPathLength.push_back((pathLenghDistInt.interpolation(M_PI/2.0)           + pathLenghDistInt.interpolation(3.0*M_PI/2.0))     /2.0);

	interpolationF temperatureDistInt(phiGridPts, temperatureDist);
	avgTemp.push_back(cubicIntegrate(phiGridPts, temperatureDist)/2.0/M_PI);
	avgTemp.push_back((temperatureDistInt.interpolation(phiGridPts.front()) + temperatureDistInt.interpolation(phiGridPts.back()))/2.0);
	avgTemp.push_back((temperatureDistInt.interpolation(M_PI/2.0)           + temperatureDistInt.interpolation(3.0*M_PI/2.0))     /2.0);
}