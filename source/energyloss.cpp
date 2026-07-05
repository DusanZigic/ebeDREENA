#include "energyloss.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <random>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "utils.hpp"
#include "polyintegrator.hpp"
#include "linearinterpolator.hpp"
#include "grids.hpp"

EnergyLoss::EnergyLoss(const config::energyLossConfig &cfg)
{
	m_collsys    = cfg.collsys;
    m_sNN        = cfg.sNN;
    m_pName      = cfg.pName;
    m_centrality = cfg.centrality;
    m_xB         = cfg.xB;
    m_BCPP       = cfg.BCPP;
    m_eventN     = cfg.eventN;
    m_phiGridN   = cfg.phiGridN;
    m_TIMESTEP   = cfg.TIMESTEP;
	m_TCRIT      = cfg.TCRIT;
    m_BCPSEED    = cfg.BCPSEED;

    m_nf = m_sNN == "200GeV" ? 2.5 : 3.0;
	double mu = utils::debyeMass(m_nf, 3.0/2.0*m_TCRIT);
	m_mgC = mu*utils::constants::INV_SQRT_2;
	if (m_pName == "Bottom") {
		m_MC = 4.75;
	} else if (m_pName == "Charm") {
		m_MC = 1.2;
	} else if (m_pName == "Gluon") {
		m_MC = mu*utils::constants::INV_SQRT_2;
	} else {
		m_MC = mu*utils::constants::INV_SQRT_6;
	}
	m_TCollConst = 3.0/2.0*m_TCRIT;
}

EnergyLoss::~EnergyLoss() {}

void EnergyLoss::runEnergyLoss()
{
	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

	if (loadLdndx()        != 0) return;
	if (loadLNorm()        != 0) return;
	if (loadLColl()        != 0) return;
	if (generateTempGrid() != 0) return;
 	if (loadPhiPoints()    != 0) return;

	if ((m_pName == "Bottom") || (m_pName == "Charm")) {
		runELossHeavyFlavour();
	}
	else if (m_pName == "LQuarks") {
		runELossLightQuarks();
	}
	else {
		runELossLightFlavour();
	}
}

int EnergyLoss::loaddsdpti2(const std::string &pname, LinearInterpolator<double> &dsdpti2int) 
{
	const std::string path_in = "./ptDists/ptDist" + m_sNN + "/ptDist_" + m_sNN + "_" + pname + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file." << std::endl;
		return 1;
	}

	std::vector<double> pTdistX, pTdistF;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; pTdistX.push_back(buffer);
		ss >> buffer; pTdistF.push_back(buffer);
	}

	dsdpti2int.setData(pTdistX, pTdistF);

	file_in.close();

	return 0;
}

int EnergyLoss::loadLdndx()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/ldndx_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open Ldndx table file." << std::endl;
		return 1;
	}

	std::vector<double> Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; Ldndx_tau.push_back(buffer);
		ss >> buffer; Ldndx_p.push_back(buffer);
		ss >> buffer; Ldndx_T.push_back(buffer);
		ss >> buffer; Ldndx_x.push_back(buffer);
		ss >> buffer; Ldndx_f.push_back(buffer);
	}

	file_in.close();

	m_Ldndx.setData(Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f);

	// TODO: domain check
	// std::vector<std::vector<double>> domain = m_Ldndx.domain();
	// if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.xPts(0)    < domain[3][0]) {std::cerr << "Error:   x grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.xPts(-1)   > domain[3][1]) {std::cerr << "Error:   x grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}

	return 0;
}

int EnergyLoss::loadLNorm()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lnorm_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LNorm table file." << std::endl;
		return 1;
	}

	std::vector<double> LNorm_tau, LNorm_p, LNorm_T, LNorm_f; //defining vectors that store LNorm table values

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; LNorm_tau.push_back(buffer);
		ss >> buffer; LNorm_p.push_back(buffer);
		ss >> buffer; LNorm_T.push_back(buffer);
		ss >> buffer; LNorm_f.push_back(buffer);
	}

	file_in.close();

	m_LNorm.setData(LNorm_tau, LNorm_p, LNorm_T, LNorm_f);

	// TODO: domain check
	// std::vector<std::vector<double>> domain = m_LNorm.domain();
	// if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}

	return 0;
}

int EnergyLoss::loadLColl()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lcoll_nf=" + nfss.str() + "_" + partName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LColl table file." << std::endl;
		return 1;
	}

	std::vector<double> LColl_p, LColl_T, LColl_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;
            
		std::stringstream ss(line);
		ss >> buffer; LColl_p.push_back(buffer);
		ss >> buffer; LColl_T.push_back(buffer);
		ss >> buffer; LColl_f.push_back(buffer);
	}

	file_in.close();

	m_LColl.setData(LColl_p, LColl_T, LColl_f);

	// TODO: domain check
	// std::vector<std::vector<double>> domain = m_LColl.domain();
	// if (m_Grids.pCollPts(0)  < domain[0][0]) {std::cerr << "Error: p grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.pCollPts(-1) > domain[0][1]) {std::cerr << "Error: p grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TCollPts(0)  < domain[1][0]) {std::cerr << "Error: T grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	// if (m_Grids.TCollPts(-1) > domain[1][1]) {std::cerr << "Error: T grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}

	return 0;
}

int EnergyLoss::generateTempGrid()
{
    const std::string path_in = "./evols/evols_cent=" + m_centrality + "/evolgridparams.dat";

    std::ifstream file_in(path_in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open evolution grid parameters file." << std::endl;
        return 1;
    }

    // NOTE (dusan): helper lambda to get the next line containing actual data
    auto readNextDataLine = [](std::ifstream &file, std::string &line) -> bool {
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '#') {
                return true;
            }
        }
        return false;
    };

    std::string line;

    // NOTE (dusan): parse and generate tau grid
    if (!readNextDataLine(file_in, line)) {
        std::cerr << "Error: failed to read tau parameters from 'evolgridparams' file." << std::endl;
        return 1;
    }
    {
        std::stringstream ss(line);
        ss >> m_tempTau0 >> m_tempTauStep;
        m_tau0 = m_tempTau0;

        m_tempTauGrid.clear();
        std::size_t i = 0;
        while (true) {
            double current_tau = m_tempTau0 + i * m_tempTauStep;
            if (current_tau >= (m_tauMaxFM + m_tempTauStep - 1e-9)) break;            
            m_tempTauGrid.push_back(current_tau);
            i++;
        }
        m_tempTauMax = m_tempTauGrid.size();
    }

    // NOTE (dusan): parse and generate x grid
    if (!readNextDataLine(file_in, line)) {
        std::cerr << "Error: failed to read x parameters from 'evolgridparams' file." << std::endl;
        return 1;
    }
    {
        double xMax;
        std::stringstream ss(line);
        ss >> m_tempX0 >> xMax >> m_tempXStep;

        m_tempXGrid.clear();
        m_tempXMax = static_cast<std::size_t>(std::round((xMax - m_tempX0) / m_tempXStep)) + 1;
        
        m_tempXGrid.reserve(m_tempXMax);
        for (std::size_t i = 0; i < m_tempXMax; ++i) {
            m_tempXGrid.push_back(m_tempX0 + i * m_tempXStep);
        }
    }

    // NOTE (dusan): parse and generate y grid
    if (!readNextDataLine(file_in, line)) {
        std::cerr << "Error: Failed to read y parameters from 'evolgridparams' file." << std::endl;
        return 1;
    }
    {
        double yMax;
        std::stringstream ss(line);
        ss >> m_tempY0 >> yMax >> m_tempYStep;

        m_tempYGrid.clear();
        m_tempYMax = static_cast<std::size_t>(std::round((yMax - m_tempY0) / m_tempYStep)) + 1;
        
        m_tempYGrid.reserve(m_tempYMax);
        for (std::size_t i = 0; i < m_tempYMax; ++i) {
            m_tempYGrid.push_back(m_tempY0 + i * m_tempYStep);
        }
    }

    file_in.close();
    return 0;
}

int EnergyLoss::loadPhiPoints()
{
	const std::string path_in = "./phiGaussPts/phiptsgauss" + std::to_string(m_phiGridN) + ".dat";
	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open phi points file. Aborting..." << std::endl;
		return 1;
	}

	std::string line; double buffer;

	while(std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> buffer; m_phiGridPts.push_back(buffer);
	}

	file_in.close();

	if (m_phiGridN != m_phiGridPts.size()) {
		std::cerr << "Error: phiGridN not equal to number of point imported from a file. Aborting..." << std::endl;
		return 1;
	}

	return 0;
}

int EnergyLoss::loadBinCollPoints(std::size_t event_id, std::vector<std::vector<double>> &bcpoints)
{
	const std::string path_in = "./binarycollpts/binarycollpts_cent=" + m_centrality + "/binarycollpts" + std::to_string(event_id) + ".dat";

	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open binary collision points file for event: " + std::to_string(event_id) + "." << std::endl;
		return 1;
	}

	std::string line; double buffer;

	std::vector<double> xpoints, ypoints;

	while (std::getline(file_in, line))
	{
		if (line.length() == 0)
            continue;
        
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; xpoints.push_back(buffer);
		ss >> buffer; ypoints.push_back(buffer);
	}

	bcpoints.resize(xpoints.size());

	for (std::size_t iBCP=0; iBCP<xpoints.size(); iBCP++)
	{
		bcpoints[iBCP].push_back(xpoints[iBCP]);
        bcpoints[iBCP].push_back(ypoints[iBCP]);
	}

	file_in.close();

	return 0;
}

int EnergyLoss::generateInitPosPoints(std::size_t event_id, std::vector<double> &xPoints, std::vector<double> &yPoints)
{
	std::vector<std::vector<double>> bcpts; if (loadBinCollPoints(event_id, bcpts) != 0) return 1;

	std::size_t bsptsNum = static_cast<std::size_t>(m_BCPP*static_cast<double>(bcpts.size()));

	if (bsptsNum < 1) bsptsNum = 1;

	if (m_BCPSEED == 0) {
		std::random_device rd; auto rng = std::default_random_engine{rd()};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}
	else {
		auto rng = std::default_random_engine{static_cast<long unsigned int>(m_BCPSEED)};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}

	for (std::size_t iBCP=0; iBCP<bsptsNum; iBCP++) {
		xPoints.push_back(bcpts[iBCP][0]);
        yPoints.push_back(bcpts[iBCP][1]);
	}

	return 0;
}

int EnergyLoss::loadTProfile(std::size_t event_id, LinearInterpolator<double> &tempProfile)
{
    const std::string path_in = "./evols/evols_cent=" + m_centrality + "/tempevol" + std::to_string(event_id) + ".dat";

    std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open temperature evolution file for event " + std::to_string(event_id) + "." << std::endl;
        return 1;
    }

    std::vector<double> temps; 
    float buffer;

    while (file_in.read(reinterpret_cast<char*>(&buffer), sizeof(buffer))) { // NOTE (dusan): reading 32bit floats from binary data file
        temps.push_back(static_cast<double>(buffer));
    }

    file_in.close();

    std::size_t spatialGridSize = m_tempXMax * m_tempYMax;
    if (temps.size() % spatialGridSize != 0) {
        std::cerr << "Error: data size is not a perfect multiple of the spatial grid layout." << std::endl;
        return 1;
    }

    std::size_t currentTauN = temps.size() / spatialGridSize;

    if (currentTauN > m_tempTauGrid.size()) {
        std::cerr << "Error: imported profile's tau length larger than maximum pre-allocated grid size." << std::endl;
        return 1;
    }

    tempProfile.setData(
        std::vector<double>(m_tempTauGrid.begin(), m_tempTauGrid.begin() + currentTauN),
        m_tempXGrid,
        m_tempYGrid,
        temps
    );

    return 0;
}

void EnergyLoss::generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab) const
//function that generates sampling points for Gaussian integration
//qGTab, fGTab - vectors that store sampling point <- output
{	
	double sigmaNum = 3.5; //setting sigma
	double sigmaStep = 0.25; //setting step
	std::size_t GTabLen = 2 * static_cast<std::size_t>(sigmaNum / sigmaStep) + 1; //setting length of sampling points
	
	double GaussTabSum = 0.0; //setting normalization sum to zero
	
	for (std::size_t iG=0; iG<GTabLen; iG++) //calculating sampling points
	{
		qGTab.push_back(-1.0*sigmaNum + static_cast<double>(iG)*sigmaStep); //setting qGaussTab values
		fGTab.push_back(std::exp(-qGTab.back()*qGTab.back()/2.0));          //setting fGaussTab values
		GaussTabSum += fGTab.back();                                        //adding to normalization sum
	}
	
	for (std::size_t iG=0; iG<GTabLen; iG++)  //normalizing
	{
		fGTab[iG] /= GaussTabSum; //dividing fGaussTab values with total sum
	}
}

void EnergyLoss::calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp) const
{
	LinearInterpolator<double> pathLenghDistInt(m_phiGridPts, pathLenghDist);
	avgPathLength.push_back(poly::cubicIntegrate(m_phiGridPts, pathLenghDist)/2.0/utils::constants::PI);
	avgPathLength.push_back((pathLenghDistInt.interpolate(m_phiGridPts.front())     + pathLenghDistInt.interpolate(m_phiGridPts.back()))          / 2.0);
	avgPathLength.push_back((pathLenghDistInt.interpolate(utils::constants::PI/2.0) + pathLenghDistInt.interpolate(3.0*utils::constants::PI/2.0)) / 2.0);

	LinearInterpolator<double> temperatureDistInt(m_phiGridPts, temperatureDist);
	avgTemp.push_back(poly::cubicIntegrate(m_phiGridPts, temperatureDist)/2.0/utils::constants::PI);
	avgTemp.push_back((temperatureDistInt.interpolate(m_phiGridPts.front())      + temperatureDistInt.interpolate(m_phiGridPts.back()))          / 2.0);
	avgTemp.push_back((temperatureDistInt.interpolate(utils::constants::PI/2.0)  + temperatureDistInt.interpolate(3.0*utils::constants::PI/2.0)) / 2.0);
}

int EnergyLoss::exportResults(const std::string &particleName, std::size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, std::size_t trajecNum, std::size_t elossNum) const
{
	std::vector<std::string> header;
    header.push_back("#collision_system: " + m_collsys);
	header.push_back("#collision_energy: " + m_sNN);
	header.push_back("#particle_type: " + particleName);
	header.push_back("#centrality: " + m_centrality);

	std::stringstream xbsstr; xbsstr << std::fixed << std::setprecision(1) << m_xB;
	header.push_back("#xB = " + xbsstr.str());

	header.push_back("#event_id: " + std::to_string(event_id));

	std::stringstream avgPathLengthSStr[3];
    for (std::size_t i=0; i<3; i++) avgPathLengthSStr[i] << std::fixed << std::setprecision(6) << avgPathLength[i];
	header.push_back("#average_path-lengths: " + avgPathLengthSStr[0].str() + ", " + avgPathLengthSStr[1].str() + ", " + avgPathLengthSStr[2].str());

	std::stringstream avgTempSStr[3];
    for (std::size_t i=0; i<3; i++) avgTempSStr[i] << std::fixed << std::setprecision(6) << avgTemp[i];
	header.push_back("#average_temperatures: " + avgTempSStr[0].str() + ", " + avgTempSStr[1].str() + ", " + avgTempSStr[2].str());
	
	header.push_back("#number_of_angles:                " + std::to_string(m_phiGridN));
	
	header.push_back("#total_number_of_trajectories:    " + std::to_string(trajecNum));
	header.push_back("#total_number_of_jet_energy_loss: " + std::to_string(elossNum));

	header.push_back("#BCPSEED: " + std::to_string(m_BCPSEED));

	header.push_back("#-------------------------------------------------------");
	header.push_back("#   pT [GeV]       phi          R_AA   ");

	//setting file path:
	const std::string path_out = "./results/results" + particleName + "/" + particleName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xbsstr.str() + "_dist_" + std::to_string(event_id) + ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable to open RAA(pT,phi) distribution file for event " + std::to_string(event_id) + "." << std::endl;
		return 1;
	}

	for (const auto &h : header) file_out << h << "\n";

	for (std::size_t ipT= 0; ipT<m_Grids.finPtsLength(); ipT++) //printing RAA(pT,phi) to file
		for (std::size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
			file_out << std::fixed << std::setw(14) << std::setprecision(10) <<   m_Grids.finPts(ipT) << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) <<    m_phiGridPts[iPhi] << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) << RAApTphi[ipT][iPhi] << "\n";
		}

	file_out.close();

	return 0;
}

void EnergyLoss::RadCollEL(double X0, double Y0, double phi0,
						   const LinearInterpolator<double> &TProfile,
						   std::vector<double> &radiativeRAA1, std::vector<std::vector<double>> &radiativeRAA2,
						   std::vector<double> &collisionalEL,
						   double &pathLength, double &temp) const noexcept
{
	std::vector<double> currLTTabL, currLTTabT; // NOTE (dusan): holds tau (L) and T for a given trajectory
	currLTTabL.reserve(256);                    // NOTE (dusan): default reservation to prevent vector growth overhead
    currLTTabT.reserve(256);

	double t = m_tau0;
    double currTemp;
    const double cos_phi0 = std::cos(phi0);
    const double sin_phi0 = std::sin(phi0);

	// NOTE (dusan): profile trajectory propagation
	while ((currTemp = TProfile.interpolate(t, X0 + t * cos_phi0, Y0 + t * sin_phi0)) > m_TCRIT) {
        currLTTabL.push_back(t);
        currLTTabT.push_back(currTemp);
        t += m_TIMESTEP;
    }

	if (currLTTabL.size() > 1) { // NOTE (dusan): calculating energy loss if path-length is longer than thermalization time
		const std::size_t numL = currLTTabL.size();
        const std::size_t numP = m_Grids.pPts().size();
        const std::size_t numX = m_Grids.xPts().size();

		std::vector<double> currNormF(numL);   // NOTE (dusan): norm = norm(L) for current trajectory and p
        std::vector<double> normSparseP(numP); // NOTE (dusan): norm = norm(p) integrated over L
		std::vector<double> normSparseF(numP);

		std::vector<double> currDndxF(numL);          // NOTE (dusan): dndx = dndx(L) for current trajectory and p, x
        std::vector<double> dndxSparseP(numP * numX); // NOTE (dusan): dndx = dndx(p, x) integrated over L
		std::vector<double> dndxSparseX(numP * numX);
		std::vector<double> dndxSparseF(numP * numX);

		for (std::size_t iP = 0; iP < numP; ++iP) {
			const double p = m_Grids.pPts()[iP];

			for (std::size_t iL = 0; iL < numL; ++iL) {
        		currNormF[iL] = m_LNorm.interpolate(currLTTabL[iL], p, currLTTabT[iL]);
    		}
			normSparseP[iP] = p;
    		normSparseF[iP] = poly::linearIntegrate(currLTTabL, currNormF);

			for (std::size_t ix = 0; ix < numX; ++ix) {
				const double x = m_Grids.xPts()[ix];
				for (std::size_t iL = 0; iL < numL; ++iL) {
        		    currDndxF[iL] = m_Ldndx.interpolate(currLTTabL[iL], p, currLTTabT[iL], x);
        		}
				dndxSparseP[iP * numX + ix] = p;
        		dndxSparseX[iP * numX + ix] = x;
        		dndxSparseF[iP * numX + ix] = poly::linearIntegrate(currLTTabL, currDndxF);
			}
		}

		LinearInterpolator<double> currNorm(normSparseP, normSparseF);
        LinearInterpolator<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseF);

		const std::size_t numRad = m_Grids.RadPts().size();
        const std::size_t numFdp = m_Grids.FdpPts().size();

		radiativeRAA1.reserve(numRad);
        radiativeRAA2.reserve(numRad);

		for (const auto &ph : m_Grids.RadPts()) {
			radiativeRAA1.push_back(dAp410(ph, currNorm));

			radiativeRAA2.emplace_back();
			radiativeRAA2.back().reserve(numFdp);
			for (const auto &Fdp : m_Grids.FdpPts()) {
        		radiativeRAA2.back().push_back(FdA(ph, Fdp, currNorm, currDndx));
    		}
		}

		std::vector<double> currCollF(numL);
        const std::size_t numPColl = m_Grids.pCollPts().size();
        collisionalEL.reserve(numPColl);

		for (const auto &p : m_Grids.pCollPts()) 
        {
            for (std::size_t iL = 0; iL < numL; ++iL) {
                currCollF[iL] = m_LColl.interpolate(p, currLTTabT[iL]);
            }
            collisionalEL.push_back(poly::linearIntegrate(currLTTabL, currCollF));
        }

		pathLength = currLTTabL.back(); // NOTE (dusan): average path-length and temperature
        double sumTemp = 0.0;
        for (const double tVal : currLTTabT) {
            sumTemp += tVal;
        }
        temp = sumTemp / static_cast<double>(numL);
	}

	else { // NOTE (dusan): path-length is smaller than thermalization time

		pathLength = 0.0;
		temp       = 0.0;
	}
}

void EnergyLoss::RadCollEL(double X0, double Y0, double phi0,
						   const LinearInterpolator<double> &TProfile,
						   std::vector<double> &radiativeRAA,
						   std::vector<double> &collisionalEL,
						   double &pathLength, double &temp) const noexcept
{
	std::vector<double> currLTTabL, currLTTabT; // NOTE (dusan): holds tau (L) and T for a given trajectory
	currLTTabL.reserve(256);                    // NOTE (dusan): default reservation to prevent vector growth overhead
    currLTTabT.reserve(256);

	double t = m_tau0;
    double currTemp;
    const double cos_phi0 = std::cos(phi0);
    const double sin_phi0 = std::sin(phi0);

	// NOTE (dusan): profile trajectory propagation
	while ((currTemp = TProfile.interpolate(t, X0 + t * cos_phi0, Y0 + t * sin_phi0)) > m_TCRIT) {
        currLTTabL.push_back(t);
        currLTTabT.push_back(currTemp);
        t += m_TIMESTEP;
    }
	
	if (currLTTabL.size() > 1) { // NOTE (dusan): calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		const std::size_t numL = currLTTabL.size();
        const std::size_t numP = m_Grids.pPts().size();
        const std::size_t numX = m_Grids.xPts().size();

		std::vector<double> currNormF(numL);   // NOTE (dusan): norm = norm(L) for current trajectory and p
        std::vector<double> normSparseP(numP); // NOTE (dusan): norm = norm(p) integrated over L
		std::vector<double> normSparseF(numP);

		std::vector<double> currDndxF(numL);          // NOTE (dusan): dndx = dndx(L) for current trajectory and p, x
        std::vector<double> dndxSparseP(numP * numX); // NOTE (dusan): dndx = dndx(p, x) integrated over L
		std::vector<double> dndxSparseX(numP * numX);
		std::vector<double> dndxSparseF(numP * numX);

		for (std::size_t iP = 0; iP < numP; ++iP) {
			const double p = m_Grids.pPts()[iP];

			for (std::size_t iL = 0; iL < numL; ++iL) {
        		currNormF[iL] = m_LNorm.interpolate(currLTTabL[iL], p, currLTTabT[iL]);
    		}
			normSparseP[iP] = p;
    		normSparseF[iP] = poly::linearIntegrate(currLTTabL, currNormF);

			for (std::size_t ix = 0; ix < numX; ++ix) {
				const double x = m_Grids.xPts()[ix];
				for (std::size_t iL = 0; iL < numL; ++iL) {
        		    currDndxF[iL] = m_Ldndx.interpolate(currLTTabL[iL], p, currLTTabT[iL], x);
        		}
				dndxSparseP[iP * numX + ix] = p;
        		dndxSparseX[iP * numX + ix] = x;
        		dndxSparseF[iP * numX + ix] = poly::linearIntegrate(currLTTabL, currDndxF);
			}
		}

		LinearInterpolator<double> currNorm(normSparseP, normSparseF);
        LinearInterpolator<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseF);
		
		radiativeRAA.reserve(m_Grids.RadPts().size());
		
		for (const auto &p : m_Grids.RadPts())
			radiativeRAA.push_back(dA41(p, currNorm, currDndx) / m_dsdpti2.interpolate(p));
		
		std::vector<double> currCollF(numL);
        const std::size_t numPColl = m_Grids.pCollPts().size();
        collisionalEL.reserve(numPColl);

		for (const auto &p : m_Grids.pCollPts()) 
        {
            for (std::size_t iL = 0; iL < numL; ++iL) {
                currCollF[iL] = m_LColl.interpolate(p, currLTTabT[iL]);
            }
            collisionalEL.push_back(poly::linearIntegrate(currLTTabL, currCollF));
        }

		pathLength = currLTTabL.back(); // NOTE (dusan): average path-length and temperature
        double sumTemp = 0.0;
        for (const double tVal : currLTTabT) {
            sumTemp += tVal;
        }
        temp = sumTemp / static_cast<double>(numL);
	}
	
	else { // NOTE (dusan): path-length is smaller than thermalization time

		pathLength = 0.0;
		temp       = 0.0;
	}
}

void EnergyLoss::runELossHeavyFlavour()
{
	if (loaddsdpti2(m_pName, m_dsdpti2) != 0) return;

	FdAHaltonSeqInit(150);

	#pragma omp parallel for schedule(dynamic)
	for (std::size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		LinearInterpolator<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<double>> RAAdist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		std::size_t trajectoryNum = 0, energylossNum = 0;

		for (std::size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<double> sumRAA1(m_Grids.finPtsLength(), 0.0);

			std::vector<std::vector<double>> sumRAA2(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0));

			std::size_t pltCNT = 0; //path-length and temperature distribution counter

			for (std::size_t iXY=0; iXY<xPoints.size(); iXY++)
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++;

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<double> singleRAA1; std::vector<std::vector<double>> singleRAA2;
					gaussFilterIntegrate(radRAA1, radRAA2, collEL, singleRAA1, singleRAA2);

					for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
						sumRAA1[iFinPts] += singleRAA1[iFinPts];
						for (std::size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++)
							sumRAA2[iFinPts][iFdp] += singleRAA2[iFinPts][iFdp];
					}
				}
				else { //if path length is smaller than tau0:

					for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) //adding RAA1, which is 1.0, to RAA sum; RAA2 is 0 in this case
						sumRAA1[iFinPts] += 1.0; 
				}
			}

			double weightsum = static_cast<double>(xPoints.size());
			std::for_each(sumRAA1.begin(), sumRAA1.end(), [weightsum](double &c){ c/=weightsum; });
			for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				std::for_each(sumRAA2[iFinPts].begin(), sumRAA2[iFinPts].end(), [weightsum](double &c){ c/=weightsum; });

			//setting RAA(pT,phi) value by integrating over p:
			for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				RAAdist[iFinPts][iPhi] = sumRAA1[iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iFinPts])/m_Grids.finPts(iFinPts);

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
		
		exportResults(m_pName, eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void EnergyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    LinearInterpolator<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        LinearInterpolator<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolate(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

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
			for (std::size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (m_dsdpti2.interpolate(pT + dppT)*RadRelInt.interpolate(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / m_dsdpti2.interpolate(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		LinearInterpolator<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolate(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();            //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (std::size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (m_dsdpti2.interpolate(pT + dpT + dppT)*RadRelInt.interpolate(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / m_dsdpti2.interpolate(pT) * GFSum);
			}
		}
	}
}

void EnergyLoss::runELossLightQuarks()
{
	const std::vector<std::string> lightQuarksList{"Down", "DownBar", "Strange", "Up", "UpBar"};

	std::vector<LinearInterpolator<double>> dsdpti2LightQuarks(lightQuarksList.size());

	for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
		if (loaddsdpti2(lightQuarksList[iLQ], dsdpti2LightQuarks[iLQ]) != 0) return;

	FdAHaltonSeqInit(100);

	#pragma omp parallel for schedule(dynamic)
	for (std::size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		LinearInterpolator<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<std::vector<double>>> RAAdist(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0)));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		std::size_t trajectoryNum = 0, energylossNum = 0;

		for (std::size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<std::vector<double>> sumRAA1(lightQuarksList.size(), std::vector<double>(m_Grids.finPtsLength(), 0.0));

			std::vector<std::vector<std::vector<double>>> sumRAA2(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0)));

			std::size_t pltCNT = 0; //path-length and temperature distribution counter

			for (std::size_t iXY=0; iXY<xPoints.size(); iXY++) //loop over x and y initial position points
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++; //adding to number of energy loss calculations

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<std::vector<double>> singleRAA1(lightQuarksList.size());
					std::vector<std::vector<std::vector<double>>> singleRAA2(lightQuarksList.size());

					for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						gaussFilterIntegrate(dsdpti2LightQuarks[iLQ], radRAA1, radRAA2, collEL, singleRAA1[iLQ], singleRAA2[iLQ]);

					for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
						for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
							sumRAA1[iLQ][iFinPts] += singleRAA1[iLQ][iFinPts];
							for (std::size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++)
								sumRAA2[iLQ][iFinPts][iFdp] += singleRAA2[iLQ][iFinPts][iFdp];
						}
					}
				}
				else {
					for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
							sumRAA1[iLQ][iFinPts] += 1.0;
				}
			}

			double weightsum = static_cast<double>(xPoints.size());
			for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
				std::for_each(sumRAA1[iLQ].begin(), sumRAA1[iLQ].end(), [weightsum](double &c){ c/=weightsum; });
				for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
					std::for_each(sumRAA2[iLQ][iFinPts].begin(), sumRAA2[iLQ][iFinPts].end(), [weightsum](double &c){ c/=weightsum; });
			}

			//setting RAA(pT,phi) value by integrating over p:
			for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
				for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
					RAAdist[iLQ][iFinPts][iPhi] = sumRAA1[iLQ][iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iLQ][iFinPts])/m_Grids.finPts(iFinPts);

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
		
		for (std::size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
			exportResults(lightQuarksList[iLQ], eventID, RAAdist[iLQ], avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void EnergyLoss::gaussFilterIntegrate(const LinearInterpolator<double> &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
//dsdpti2lquark - light quark initial pT distribution      						  <- input
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    LinearInterpolator<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        LinearInterpolator<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolate(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

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
			for (std::size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (dsdpti2lquark.interpolate(pT + dppT)*RadRelInt.interpolate(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / dsdpti2lquark.interpolate(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		LinearInterpolator<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolate(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (std::size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (dsdpti2lquark.interpolate(pT + dpT + dppT)*RadRelInt.interpolate(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2lquark.interpolate(pT) * GFSum);
			}
		}
	}
}

void EnergyLoss::runELossLightFlavour()
{
	if (loaddsdpti2(m_pName, m_dsdpti2) != 0) return;

	dAHaltonSeqInit(1000);

	#pragma omp parallel for schedule(dynamic)
	for (std::size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		LinearInterpolator<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<double>> RAAdist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		std::size_t trajectoryNum = 0, energylossNum = 0;

		for (std::size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<double> sumRAA(m_Grids.finPtsLength(), 0.0);

			std::size_t pltCNT = 0; //path-length and temperature distribution counter

			for (std::size_t iXY=0; iXY<xPoints.size(); iXY++)
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA, collEL; double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++; //adding to number of energy loss calculations

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<double> singleRAA;
					gaussFilterIntegrate(radRAA, collEL, singleRAA);

					for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
						sumRAA[iFinPts] += singleRAA[iFinPts];
				}
				else { //if path length is smaller than tau0:

					for (std::size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
						sumRAA[iFinPts] += 1.0;
				}
			}

			double weightsum = (double)(xPoints.size());
			for (std::size_t iFinPts= 0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				RAAdist[iFinPts][iPhi] = sumRAA[iFinPts]/weightsum;

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);

		exportResults(m_pName, eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void EnergyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA) const
//function that performs Gauss filter integration - default algorithm
//radiativeRAA  - raditive RAA 							   <- input
//collisionalEL - collisional energy loss				   <- input
//singRAA 		- RAA array after Gauss filter integration <- output
{
    LinearInterpolator<double> RadRelInt(m_Grids.RadPts(),   radiativeRAA);  //creating radiative RAA interpolated function
    LinearInterpolator<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	double GFSum; //defining sum variable for Gauss filter

	double dpT; //defining pT and dpT variables

	double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value

	double sigmaColl; //defining variable for collisional sigma
	
	//Gauss filter
	for (const auto &pT : m_Grids.finPts())
	{
		GFSum = 0.0L;

		muCollCurrVal = muCollInt.interpolate(pT);

		sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

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
		for (std::size_t iG=0; iG<qGaussTab.size(); iG++)
		{
			dpT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
			GFSum += (m_dsdpti2.interpolate(pT + dpT)*RadRelInt.interpolate(pT + dpT)*(pT + dpT) / pT * fGaussTab[iG]);
		}

		singRAA.push_back(1.0 / m_dsdpti2.interpolate(pT) * GFSum);
	}
}