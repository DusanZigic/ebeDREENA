#include "importexport.hpp"
#include "arsenal.hpp"
#include "grids.hpp"
#include "linearinterpolation.hpp"
#include "elossheader.hpp"
// #include "ltables.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <cctype>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cctype>

int loaddsdpti2()
{
	const std::string path_in = "./ptDists/ptDist" + sNN + "/ptDist_" + sNN + "_" + pName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file." << std::endl;
		return -1;
	}

	std::vector<double> pTdistX, pTdistF; //defining vectors that store dsdpti2 values

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; pTdistX.push_back(buffer);
		ss >> buffer; pTdistF.push_back(buffer);
	}

	dsdpti2.SetData(pTdistX, pTdistF);

	file_in.close();

	return 1;
}

int loaddsdpti2(const std::string &pname, interpolationF &dsdpti2int)
{
	const std::string path_in = "./ptDists/ptDist" + sNN + "/ptDist_" + sNN + "_" + pname + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file." << std::endl;
		return -1;
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

	dsdpti2int.SetData(pTdistX, pTdistF);

	file_in.close();

	return 1;
}

int loadLdndx()
{
	std::string partName;
	if (pName == "Bottom") partName = "Bottom";
	else if (pName == "Charm") partName = "Charm";
	else if (pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << nf;

	const std::string path_in = "./LTables/LdndxTbl_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open Ldndx table file." << std::endl;
		return -1;
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

	Ldndx.SetData(Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f);

	file_in.close();

	return 1;
}

int loadLNorm()
{
	std::string partName;
	if (pName == "Bottom") partName = "Bottom";
	else if (pName == "Charm") partName = "Charm";
	else if (pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << nf;

	const std::string path_in = "./LTables/LNormTbl_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LNorm table file." << std::endl;
		return -1;
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

	LNorm.SetData(LNorm_tau, LNorm_p, LNorm_T, LNorm_f);

	file_in.close();

	return 1;
}

int loadLColl()
{
	std::string partName;
	if (pName == "Bottom") partName = "Bottom";
	else if (pName == "Charm") partName = "Charm";
	else if (pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << nf;

	const std::string path_in = "./LTables/LCollTbl_nf=" + nfss.str() + "_" + partName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LColl table file." << std::endl;
		return -1;
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

	LColl.SetData(LColl_p, LColl_T, LColl_f);

	file_in.close();

	return 1;
}

static double tau_max = 25.0;					    //maximal value of tau for TProfile grids in fm
double temp_tau0, temp_tauStep; size_t temp_tauMax; //temperature evolution tau grid parameters
double   temp_x0,   temp_xStep; size_t   temp_xMax; //temperature evolution   x grid parameters
double   temp_y0,   temp_yStep; size_t   temp_yMax; //temperature evolution   y grid parameters

std::vector<double> tempTauGrid, tempXGrid, tempYGrid; //temperature evolutions grids

int loadTempGridParams()
{
	const std::string path_in = "./TProfiles/TProfiles_bin_cent=" + centrality + "/temp_grids.dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open TProfile grid parameters file." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer; temp_tau0    = buffer;
							        ss >> buffer; temp_tauStep = buffer;
							   			          temp_tauMax  = 0;
		double tau = temp_tau0; while (tau < (tau_max+temp_tauStep)) {temp_tauMax++; tau+=temp_tauStep;}
		tau0 = temp_tau0;
	}

	{
		double x_max;
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer;    temp_x0 = buffer;
							        ss >> buffer;      x_max = buffer;
							        ss >> buffer; temp_xStep = buffer;
		double x = temp_x0; temp_xMax = 0; while (x <= x_max) {temp_xMax++; x+=temp_xStep;}
	}

	{
		double y_max;
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer;    temp_y0 = buffer;
							        ss >> buffer;      y_max = buffer;
							        ss >> buffer; temp_yStep = buffer;
		double y = temp_x0; temp_yMax = 0; while (y <= y_max) {temp_yMax++; y+=temp_xStep;}
	}

	return 1;
}

int loadTProfile(size_t event_id, interpolationF &tempProfile)
{
	const std::string path_in = "./TProfiles/TProfiles_bin_cent=" + centrality + "/TProfile_" + std::to_string(event_id) + ".dat";

	std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open temperature evolution file for event " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

    std::vector<double> temps; float buffer;

    while (true) {
        file_in.read((char*)&buffer, sizeof(buffer));
        if (file_in.eof()) break;
        temps.push_back(static_cast<double>(buffer));
    }

	file_in.close();

	if (temps.size() > tempTauGrid.size()) {
		std::cerr << "Error: imported profile's size larger than grid size." << std::endl;
		return -2;
	}

	tempProfile.SetData(tempTauGrid, tempXGrid, tempYGrid, temps);

	return 1;
}

int loadBinCollPoints(size_t event_id, std::vector<std::vector<double>> &bcpoints)
{
	const std::string path_in = "./BinaryCollPoints/BinaryCollPoints_cent=" + centrality + "/BinaryCollPoints_" + std::to_string(event_id) + ".dat";

	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open binary collision points file for event: " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	std::vector<double> xpoints, ypoints;

	while (std::getline(file_in, line))
	{
		if (line.length() > 0)
            continue;
        
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; xpoints.push_back(buffer);
		ss >> buffer; ypoints.push_back(buffer);
	}

	bcpoints.resize(xpoints.size());

	for (size_t iBCP=0; iBCP<xpoints.size(); iBCP++)
	{
		bcpoints[iBCP].push_back(xpoints[iBCP]);
        bcpoints[iBCP].push_back(ypoints[iBCP]);
	}

	file_in.close();

	return 1;
}

int loadPhiPoints()
{
	const std::string path_in = "./phiGaussPts/phiptsgauss" + std::to_string(phiGridN) + ".dat";
	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open phi points file. Aborting..." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	while(std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> buffer; phiGridPts.push_back(buffer);
	}

	phiGridN = phiGridPts.size();

	file_in.close();

	return 1;
}

int exportResults(size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, size_t trajecNum, size_t elossNum)
{
	std::vector<std::string> header;
	header.push_back("#collision_system: Pb+Pb");
	if (sNN ==  "200GeV") header[0] = "#collision_system: Au+Au";
	if (sNN == "5440GeV") header[0] = "#collision_system: Xe+Xe";
	header.push_back("#collision_energy: " + sNN);
	header.push_back("#particle_type: " + pName);
	header.push_back("#centrality: " + centrality);

	std::stringstream xbsstr; xbsstr << std::fixed << std::setprecision(1) << xB;
	header.push_back("#xB = " + xbsstr.str());

	header.push_back("#event_id: " + std::to_string(event_id));

	std::stringstream avgPathLengthSStr[3];
    for (size_t i=0; i<3; i++) avgPathLengthSStr[i] << std::fixed << std::setprecision(6) << avgPathLength[i];
	header.push_back("#average_path-lengths: " + avgPathLengthSStr[0].str() + ", " + avgPathLengthSStr[1].str() + ", " + avgPathLengthSStr[2].str());

	std::stringstream avgTempSStr[3];
    for (size_t i=0; i<3; i++) avgTempSStr[i] << std::fixed << std::setprecision(6) << avgTemp[i];
	header.push_back("#average_temperatures: " + avgTempSStr[0].str() + ", " + avgTempSStr[1].str() + ", " + avgTempSStr[2].str());
	
	header.push_back("#number_of_angles:                " + std::to_string(phiGridN));
	
	header.push_back("#total_number_of_trajectories:    " + std::to_string(trajecNum));
	header.push_back("#total_number_of_jet_energy_loss: " + std::to_string(elossNum));

	header.push_back("#-------------------------------------------------------");
	header.push_back("#   pT [GeV]       phi          R_AA   ");

	//setting file path:
	const std::string path_out = "./CResults/CResults_" + pName + "/" + pName + "_sNN=" + sNN + "_cent=" + centrality + "_xB=" + xbsstr.str() + "_dist_" + std::to_string(event_id) + ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable to open RAA(pT,phi) distribution file for event " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

	for (const auto &head : header) file_out << head << "\n";


	for (size_t ipT= 0; ipT<Grids.finPtsLength(); ipT++)
		for (size_t iPhi=0; iPhi<phiGridN; iPhi++)
			file_out << std::fixed << std::setw(14) << std::setprecision(10) <<   Grids.finPts(ipT) << " "
					 << std::fixed << std::setw(12) << std::setprecision(10) <<    phiGridPts[iPhi] << " "
					 << std::fixed << std::setw(12) << std::setprecision(10) << RAApTphi[ipT][iPhi] << "\n";

	file_out.close();

	return 1;
}

int exportResults(const std::string &particleName, size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, size_t trajecNum, size_t elossNum)
{
		std::vector<std::string> header;
	header.push_back("#collision_system: Pb+Pb");
	if (sNN ==  "200GeV") header[0] = "#collision_system: Au+Au";
	if (sNN == "5440GeV") header[0] = "#collision_system: Xe+Xe";
	header.push_back("#collision_energy: " + sNN);
	header.push_back("#particle_type: " + particleName);
	header.push_back("#centrality: " + centrality);

	std::stringstream xbsstr; xbsstr << std::fixed << std::setprecision(1) << xB;
	header.push_back("#xB = " + xbsstr.str());

	header.push_back("#event_id: " + std::to_string(event_id));

	std::stringstream avgPathLengthSStr[3];
    for (size_t i=0; i<3; i++) avgPathLengthSStr[i] << std::fixed << std::setprecision(6) << avgPathLength[i];
	header.push_back("#average_path-lengths: " + avgPathLengthSStr[0].str() + ", " + avgPathLengthSStr[1].str() + ", " + avgPathLengthSStr[2].str());

	std::stringstream avgTempSStr[3];
    for (size_t i=0; i<3; i++) avgTempSStr[i] << std::fixed << std::setprecision(6) << avgTemp[i];
	header.push_back("#average_temperatures: " + avgTempSStr[0].str() + ", " + avgTempSStr[1].str() + ", " + avgTempSStr[2].str());
	
	header.push_back("#number_of_angles:                " + std::to_string(phiGridN));
	
	header.push_back("#total_number_of_trajectories:    " + std::to_string(trajecNum));
	header.push_back("#total_number_of_jet_energy_loss: " + std::to_string(elossNum));

	header.push_back("#-------------------------------------------------------");
	header.push_back("#   pT [GeV]       phi          R_AA   ");

	//setting file path:
	const std::string path_out = "./CResults/CResults_" + pName + "/" + pName + "_sNN=" + sNN + "_cent=" + centrality + "_xB=" + xbsstr.str() + "_dist_" + std::to_string(event_id) + ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable to open RAA(pT,phi) distribution file for event " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

	for (const auto &head : header) file_out << head << "\n";


	for (size_t ipT= 0; ipT<Grids.finPtsLength(); ipT++) //printing RAA(pT,phi) to file
		for (size_t iPhi=0; iPhi<phiGridN; iPhi++)
			file_out << std::fixed << std::setw(14) << std::setprecision(10) <<   Grids.finPts(ipT) << " "
					 << std::fixed << std::setw(12) << std::setprecision(10) <<    phiGridPts[iPhi] << " "
					 << std::fixed << std::setw(12) << std::setprecision(10) << RAApTphi[ipT][iPhi] << "\n";

	file_out.close();

	return 1;
}

/*int exportLTables(const std::vector<double> &ldndxTable, const std::vector<double> &lnormTable, const std::vector<double> &lcollTable)
{
	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << LT_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << LT_nf;

	const std::string path_out_ldndx = "./LTables/LdndxTbl_nf=" + nfss.str() + "_" + LT_pName + "_xB=" + xBss.str() + ".dat";
	std::ofstream file_out_ldndx(path_out_ldndx);
	if (!file_out_ldndx.is_open()) {
		cerr << "Error: unable to open Ldndx table export file." << endl;
		return -1;
	}

	//opening LNorm export file:
	string path_out_lnorm = "LTables/LNormTbl_nf=" + nfss.str() + "_" + LT_pName + "_xB=" + xBss.str() + ".dat";
	ofstream file_out_lnorm(path_out_lnorm);
	if (!file_out_lnorm.is_open()) {
		cerr << "Error: unable to open LNorm table export file." << endl;
		return -1;
	}

	for (int tau_i=0; tau_i<LT_Grids.tauPtsLength(); tau_i++)
	{
		for (int p_i=0; p_i<LT_Grids.pPtsLength(); p_i++)
		{
			for (int T_i=0; T_i<LT_Grids.TPtsLength(); T_i++)
			{
				for (int x_i=0; x_i<LT_Grids.xPtsLength(); x_i++)
				{
					//calculating ldndx index:
					int ldndx_index = tau_i*LT_Grids.xPtsLength()*LT_Grids.TPtsLength()*LT_Grids.pPtsLength()
									+ p_i*LT_Grids.xPtsLength()*LT_Grids.TPtsLength()
									+ T_i*LT_Grids.xPtsLength()
									+ x_i;

					//printing ldndx to file:
					file_out_ldndx << fixed << setprecision(10) << LT_Grids.tauPts(tau_i) << " "
								   << fixed << setprecision(10) << LT_Grids.pPts(p_i) << " "
								   << fixed << setprecision(10) << LT_Grids.TPts(T_i) << " "
								   << fixed << setprecision(10) << LT_Grids.xPts(x_i) << " "
								   << fixed << setprecision(10) << ldndxtbl[ldndx_index] << endl;
				}

				//calculating lnorm index:
				int lnorm_index = tau_i*LT_Grids.TPtsLength()*LT_Grids.pPtsLength()
								+ p_i*LT_Grids.TPtsLength()
								+ T_i;

				//printing lnorm to file:
				file_out_lnorm << fixed << setprecision(10) << LT_Grids.tauPts(tau_i) << " "
							   << fixed << setprecision(10) << LT_Grids.pPts(p_i) << " "
							   << fixed << setprecision(10) << LT_Grids.TPts(T_i) << " "
							   << fixed << setprecision(10) << lnormtbl[lnorm_index] << endl;
			}
		}
	}

	file_out_ldndx.close(); file_out_lnorm.close(); //closing Ldndx and LNorm files

	//opening LColl export file:
	string path_out_lcoll = "LTables/LCollTbl_nf=" + nfss.str() + "_" + LT_pName + ".dat";
	ofstream file_out_lcoll(path_out_lcoll);
	if (!file_out_lcoll.is_open()) {
		cerr << "Error: unable to open LColl table export file." << endl;
		return -1;
	}

	//printing lcoll to file:
	for (int p_i=0; p_i<LT_Grids.pCollPtsLength(); p_i++)
	{
		for (int T_i=0; T_i<LT_Grids.TCollPtsLength(); T_i++)
		{
			file_out_lcoll << fixed << setprecision(10) << LT_Grids.pCollPts(p_i) << " "
						   << fixed << setprecision(10) << LT_Grids.TCollPts(T_i) << " "
						   << fixed << setprecision(10) << lcolltbl[p_i*LT_Grids.TCollPtsLength() + T_i] << endl;
		}
	}

	file_out_lcoll.close(); //closing LColl file
}*/