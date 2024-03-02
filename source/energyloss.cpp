#include "mainheader.hpp"
#include "elossheader.hpp"
#include "arsenal.hpp"
#include "importexport.hpp"
#include "grids.hpp"
#include "linearinterpolation.hpp"
#include "polyintegration.hpp"
#include "daintegrals.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

std::string collsys;				//collision system
std::string sNN; 				    //collision energy
double nf;					        //effective number of flavours
std::string pName; 				    //particle name
std::string centrality;			    //centrality class
double xB;					        //xB value
double BCPP;				        //binary collision points percentage
size_t eventN;					    //number of events 
gridPoints Grids;			        //grid points
interpolationF dsdpti2; 			//initial pT distribution
interpolationF LNorm, Ldndx, LColl; //interpolated L tables
double tau0;				        //thermalization time
double mgC, MC;				        //constant particle and gluon masses used for dA integrals
double TCollConst;			        //constant temperature used for Gauss filter integration
size_t phiGridN;				    //phi points number
std::vector<double> phiGridPts;     //defining vectors that store initial position points and angles
double TIMESTEP, TCRIT; 	        //defining time step and critical temperature

static const double lambda = 0.2;

static double FdA(double ph, double dp, interpolationF &currNormInt, interpolationF &currDndxInt) {
	return (FdA411(ph, dp, currNormInt, currDndxInt) + FdA412(ph, dp, currNormInt, currDndxInt) + FdA413(ph, dp, currNormInt, currDndxInt) +
				FdA414(ph, dp, currNormInt, currDndxInt) + FdA415(ph, dp, currNormInt, currDndxInt));
}

static double dA41(double ph, interpolationF &currNormInt, interpolationF &currDndxInt)
{
	if (pName == "Gluon") { //gluon needs 7 dA integrals
		return (dA410(ph, currNormInt) + dA411(ph, currNormInt, currDndxInt) + dA412(ph, currNormInt, currDndxInt) +dA413(ph, currNormInt, currDndxInt) +
					dA414(ph, currNormInt, currDndxInt) + dA415(ph, currNormInt, currDndxInt) + dA416(ph, currNormInt, currDndxInt) +
					dA417(ph, currNormInt, currDndxInt));
	}
	else { //light quarks need 5 dA integrals
		return (dA410(ph, currNormInt) + dA411(ph, currNormInt, currDndxInt) + dA412(ph, currNormInt, currDndxInt) + dA413(ph, currNormInt, currDndxInt) +
					dA414(ph, currNormInt, currDndxInt) + dA415(ph, currNormInt, currDndxInt));
	}
}

static void RadCollEL(double X0, double Y0, double phi0, const interpolationF &TProfile, std::vector<double> &radiativeRAA1, std::vector<std::vector<double>> &radiativeRAA2, std::vector<double> &collisionalEL, double &pathLength, double &temp)
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (modefied pT integration algorithm)
//X0, Y0, phi0  - inital position and angle 					  		     <- input
//radiativeRAA1 - radiative RAA for single trajectory (dA410)	  		     <- output
//radiativeRAA2 - radiative RAA for single trajectory (rest of dA integrals) <- output
//collisionalEL - collisional energy loss for single trajectory   		     <- output
//pathL, temp - path-length and temperature for single trajectory 		     <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = tau0, currTemp; //defining current path-length (time) and temperature

	while ((currTemp = TProfile.interpolation(t, X0 + t*std::cos(phi0), Y0 + t*std::sin(phi0))) > TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											  //table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	  //table for currDndxInterp

		for (const auto &p : Grids.pPts()) //loop over ppts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currNormTabTau[iL] = currLTTabL[iL]; 								         //setting path-lengths
				currNormTabVal[iL] = LNorm.interpolation(currLTTabL[iL], p, currLTTabT[iL]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : Grids.xPts()) //loop over xpts
			{
				for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
				{
					currDndxTabTau[iL] = currLTTabL[iL]; 									        //setting path-lengths
					currDndxTabVal[iL] = Ldndx.interpolation(currLTTabL[iL], p, currLTTabT[iL], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF currNorm(NormSparseP, NormSparseV); 			    //constructing interpolated current norm
		interpolationF currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (const auto &ph : Grids.RadPts()) //loop over Radpts
		{
			radiativeRAA1.push_back(dAp410(ph, currNorm)); //calculating radiative energy loss for dA410

			radiativeRAA2.push_back(std::vector<double>()); //resizing radiativeRAA2 2d vector

			for (const auto &Fdp : Grids.FdpPts())
				radiativeRAA2.back().push_back(FdA(ph, Fdp, currNorm, currDndx)); //calculating radiative energy loss for rest of the dA integrals

		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto& p : Grids.pCollPts()) //loop over pCollPts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currCollTabTau[iL] = currLTTabL[iL]; 				         //setting path-lengths
				currCollTabVal[iL] = LColl.interpolation(p, currLTTabT[iL]); //setting LColl values
			}

			collisionalEL.push_back(linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLength = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (size_t iL=0; iL<currLTTabL.size(); iL++) temp += currLTTabT[iL];
		temp /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLength = 0.0; //setting path-length and temperature
		temp       = 0.0;
	}
}

static void RadCollEL(double X0, double Y0, double phi0, const interpolationF &TProfile, std::vector<double> &radiativeRAA, std::vector<double> &collisionalEL, double &pathLenght, double &temp)
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (standard algorithm)
//X0, Y0, phi0  - inital position and angle 					  <- input
//radiativeRAA  - radiative RAA for single trajectory 			  <- output
//collisionalEL - collisional energy loss for single trajectory   <- output
//pathL, temp - path-length and temperature for single trajectory <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = tau0, currTemp; //defining current path-length (time) and temperature

	while ((currTemp = TProfile.interpolation(t, X0 + t*cos(phi0), Y0 + t*sin(phi0))) > TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											 //table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	 //table for currDndxInterp

		for (const auto &p : Grids.pPts()) //loop over ppts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currNormTabTau[iL] = currLTTabL[iL]; 								         //setting path-lengths
				currNormTabVal[iL] = LNorm.interpolation(currLTTabL[iL], p, currLTTabT[iL]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : Grids.xPts()) //loop over xpts
			{
				for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
				{
					currDndxTabTau[iL] = currLTTabL[iL]; 									        //setting path-lengths
					currDndxTabVal[iL] = Ldndx.interpolation(currLTTabL[iL], p, currLTTabT[iL], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF currNorm(NormSparseP, NormSparseV); 			   //constructing interpolated current norm
		interpolationF currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (const auto &p : Grids.RadPts())
			radiativeRAA.push_back(dA41(p, currNorm, currDndx)/dsdpti2.interpolation(p)); //calculating radiative RAA
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto &p : Grids.pCollPts()) //loop over pCollPts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currCollTabTau[iL] = currLTTabL[iL]; 				         //setting path-lengths
				currCollTabVal[iL] = LColl.interpolation(p, currLTTabT[iL]); //setting LColl values
			}

			collisionalEL.push_back(linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLenght = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (size_t iL=0; iL<currLTTabL.size(); iL++) temp += currLTTabT[iL];
		temp /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLenght = 0.0; //setting path-length and temperature
		     temp  = 0.0;
	}
}

static void SetELParameters()
{
	double T = 3.0 / 2.0*TCRIT;
	double mu = 0.197*sqrt((-8.0*(6.0+nf)*M_PI*M_PI*T*T)/(2.0*nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+nf)*M_PI*M_PI*T*T)/(2.0*nf-33.0)/lambda/lambda));
	mgC = mu / std::sqrt(2.0);
	if (pName == "Bottom") MC = 4.75;
	else if (pName == "Charm") MC = 1.2;
	else if (pName == "Gluon") MC = mu/std::sqrt(2.0);
	else MC = mu/sqrt(6.0);
	TCollConst = T;
	nf = 3.0; if (sNN == "200GeV") nf = 2.5;
}

//function that calculates averaged energy loss:
void AverageEL()
{
	SetELParameters();

	Grids.setGridPoints(sNN, pName, TCRIT);

	if (loadLdndx() != 1) return;
	if (loadLNorm() != 1) return;
	if (loadLColl() != 1) return;

	if (generateTempGrid() != 1) return;
	
 	if (loadPhiPoints() != 1) return;

 	
 	if ((pName == "Bottom") || (pName == "Charm")) {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//HEAVY FLAVOR CALCULATION:

		if (loaddsdpti2()  != 1) return;

		FdAHaltonSeqInit(150);

		#pragma omp parallel for schedule(dynamic)
		for (size_t eventID=1; eventID<=eventN; eventID++)
		{
			std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

			interpolationF tProfile; loadTProfile(eventID, tProfile);

			std::vector<std::vector<double>> RAAdist(Grids.finPtsLength(), std::vector<double>(phiGridN, 0.0));

			std::vector<double> pathLenghDist(phiGridN, 0.0), temperatureDist(phiGridN, 0.0);

			size_t trajectoryNum = 0, energylossNum = 0;

			for (size_t iPhi=0; iPhi<phiGridN; iPhi++)
			{
				double phi = phiGridPts[iPhi];

				std::vector<double> sumRAA1(Grids.finPtsLength(), 0.0);

				std::vector<std::vector<double>> sumRAA2(Grids.finPtsLength(), std::vector<double>(Grids.FdpPtsLength(), 0.0));

				size_t pltCNT = 0; //path-length and temperature distribution counter

				for (size_t iXY=0; iXY<xPoints.size(); iXY++)
				{
					trajectoryNum++;

					double x = xPoints[iXY], y = yPoints[iXY];

					std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
					double pathLength, temperature;
					RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

					if (pathLength > tau0) { //checking if path-length is larger than thermalization time

						energylossNum++;

						pltCNT++;
						pathLenghDist[iPhi] += pathLength;
						temperatureDist[iPhi] += temperature;

                        for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

						std::vector<double> singleRAA1; std::vector<std::vector<double>> singleRAA2;
						gaussFilterIntegrate(radRAA1, radRAA2, collEL, singleRAA1, singleRAA2);

						for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++) {
							sumRAA1[iFinPts] += singleRAA1[iFinPts];
							for (size_t iFdp=0; iFdp<Grids.FdpPtsLength(); iFdp++)
								sumRAA2[iFinPts][iFdp] += singleRAA2[iFinPts][iFdp];
						}
					}
					else { //if path length is smaller than tau0:

						for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++) //adding RAA1, which is 1.0, to RAA sum; RAA2 is 0 in this case
                            sumRAA1[iFinPts] += 1.0; 
					}
				}

				double weightsum = static_cast<double>(xPoints.size());
				std::for_each(sumRAA1.begin(), sumRAA1.end(), [weightsum](double &c){ c/=weightsum; });
				for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
					std::for_each(sumRAA2[iFinPts].begin(), sumRAA2[iFinPts].end(), [weightsum](double &c){ c/=weightsum; });

				//setting RAA(pT,phi) value by integrating over p:
				for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
					RAAdist[iFinPts][iPhi] = sumRAA1[iFinPts] + cubicIntegrate(Grids.FdpPts(), sumRAA2[iFinPts])/Grids.finPts(iFinPts);

				pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
			}

			std::vector<double> avgPathLength, avgTemp;
			calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
			
			exportResults(eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
		}
	
	}
	else if (pName == "LQuarks") {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ALL LIGHT QUARKS CALCULATION:

		const std::vector<std::string> lightQuarksList{"Down", "DownBar", "Strange", "Up", "UpBar"};

		std::vector<interpolationF> dsdpti2LightQuarks(lightQuarksList.size());

		for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
            if (loaddsdpti2(lightQuarksList[iLQ], dsdpti2LightQuarks[iLQ]) != 1) return;

		FdAHaltonSeqInit(100);

		#pragma omp parallel for schedule(dynamic)
		for (size_t eventID=1; eventID<=eventN; eventID++)
		{
			std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

			interpolationF tProfile; loadTProfile(eventID, tProfile);

			std::vector<std::vector<std::vector<double>>> RAAdist(lightQuarksList.size(), std::vector<std::vector<double>>(Grids.finPtsLength(), std::vector<double>(phiGridN, 0.0)));

            std::vector<double> pathLenghDist(phiGridN, 0.0), temperatureDist(phiGridN, 0.0);

			size_t trajectoryNum = 0, energylossNum = 0;

			for (size_t iPhi=0; iPhi<phiGridN; iPhi++)
			{
				double phi = phiGridPts[iPhi];

				std::vector<std::vector<double>> sumRAA1(lightQuarksList.size(), std::vector<double>(Grids.finPtsLength(), 0.0));

				std::vector<std::vector<std::vector<double>>> sumRAA2(lightQuarksList.size(), std::vector<std::vector<double>>(Grids.finPtsLength(), std::vector<double>(Grids.FdpPtsLength(), 0.0)));

				size_t pltCNT = 0; //path-length and temperature distribution counter

				for (size_t iXY=0; iXY<xPoints.size(); iXY++) //loop over x and y initial position points
				{
					trajectoryNum++;

					double x = xPoints[iXY], y = yPoints[iXY];

					std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
					double pathLength, temperature;
					RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

					if (pathLength > tau0) { //checking if path-length is larger than thermalization time

						energylossNum++; //adding to number of energy loss calculations

						pltCNT++;
						pathLenghDist[iPhi] += pathLength;
						temperatureDist[iPhi] += temperature;

                        for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

						std::vector<std::vector<double>> singleRAA1(lightQuarksList.size());
						std::vector<std::vector<std::vector<double>>> singleRAA2(lightQuarksList.size());

						for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
							gaussFilterIntegrate(dsdpti2LightQuarks[iLQ], radRAA1, radRAA2, collEL, singleRAA1[iLQ], singleRAA2[iLQ]);

						for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
							for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++) {
								sumRAA1[iLQ][iFinPts] += singleRAA1[iLQ][iFinPts];
								for (size_t iFdp=0; iFdp<Grids.FdpPtsLength(); iFdp++)
									sumRAA2[iLQ][iFinPts][iFdp] += singleRAA2[iLQ][iFinPts][iFdp];
							}
						}
					}
					else {
						for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
							for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
								sumRAA1[iLQ][iFinPts] += 1.0;
					}
				}

				double weightsum = static_cast<double>(xPoints.size());
				for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
					std::for_each(sumRAA1[iLQ].begin(), sumRAA1[iLQ].end(), [weightsum](double &c){ c/=weightsum; });
					for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
						std::for_each(sumRAA2[iLQ][iFinPts].begin(), sumRAA2[iLQ][iFinPts].end(), [weightsum](double &c){ c/=weightsum; });
				}

				//setting RAA(pT,phi) value by integrating over p:
				for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
					for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
						RAAdist[iLQ][iFinPts][iPhi] = sumRAA1[iLQ][iFinPts] + cubicIntegrate(Grids.FdpPts(), sumRAA2[iLQ][iFinPts])/Grids.finPts(iFinPts);

				pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
			}

            std::vector<double> avgPathLength, avgTemp;
			calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
			
			for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
				exportResults(lightQuarksList[iLQ], eventID, RAAdist[iLQ], avgPathLength, avgTemp, trajectoryNum, energylossNum);
		}
	}
	else {
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//LIGHT FLAVOR CALCULATION (INDIVIDUAL PARTICLES):

		if (loaddsdpti2()  != 1) return;

		dAHaltonSeqInit(1000);

		#pragma omp parallel for schedule(dynamic)
		for (size_t eventID=1; eventID<=eventN; eventID++)
		{
			std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

			interpolationF tProfile; loadTProfile(eventID, tProfile);

			std::vector<std::vector<double>> RAAdist(Grids.finPtsLength(), std::vector<double>(phiGridN, 0.0));

			std::vector<double> pathLenghDist(phiGridN, 0.0), temperatureDist(phiGridN, 0.0);

			size_t trajectoryNum = 0, energylossNum = 0;

			for (size_t iPhi=0; iPhi<phiGridN; iPhi++)
			{
				double phi = phiGridPts[iPhi];

				std::vector<double> sumRAA(Grids.finPtsLength(), 0.0);

				size_t pltCNT = 0; //path-length and temperature distribution counter

				for (size_t iXY=0; iXY<xPoints.size(); iXY++)
				{
					trajectoryNum++;

					double x = xPoints[iXY], y = yPoints[iXY];

					std::vector<double> radRAA, collEL; double pathLength, temperature;
					RadCollEL(x, y, phi, tProfile, radRAA, collEL, pathLength, temperature);

					if (pathLength > tau0) { //checking if path-length is larger than thermalization time

						energylossNum++; //adding to number of energy loss calculations

						pltCNT++;
						pathLenghDist[iPhi] += pathLength;
						temperatureDist[iPhi] += temperature;

                        for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

						std::vector<double> singleRAA;
						gaussFilterIntegrate(radRAA, collEL, singleRAA);

						for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
                            sumRAA[iFinPts] += singleRAA[iFinPts];
					}
					else { //if path length is smaller than tau0:

						for (size_t iFinPts=0; iFinPts<Grids.finPtsLength(); iFinPts++)
                            sumRAA[iFinPts] += 1.0;
					}
				}

				double weightsum = (double)(xPoints.size());
				for (size_t iFinPts= 0; iFinPts<Grids.finPtsLength(); iFinPts++)
                    RAAdist[iFinPts][iPhi] = sumRAA[iFinPts]/weightsum;

				pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
			}

			std::vector<double> avgPathLength, avgTemp;
			calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);

			exportResults(eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
		}
	}
}