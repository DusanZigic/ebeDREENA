#include "inputparser.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>

int GetInputs(int argc, char const *argv[], std::string &csys, std::string &snn, std::string &pname, std::string &cent, double &xb, size_t &eventn, double &bcpp, size_t &phiptsn, double &tstep, double &tcrit, int &bcpseed)
{
	std::vector<std::string> inputs; for (int i=2; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: --collsys=PbPb --sNN=5020GeV --pName=Charm --centrality=30-40% --xB=0.6 --eventN=1000 --BCPP=20% --phiGridN=25 --TIMESTEP=0.1 --TCRIT=0.155 --BCPSEED=0" << std::endl;
		return 0;
	}

	std::map<std::string, std::string> inputparams;
	for (const auto &in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

	//checking if configuration file is provided:
	std::map<std::string, std::string> inputparams_f;
	if (inputparams.count("c") > 0) {
		std::ifstream file_in(inputparams["c"]);
		if (!file_in.is_open()) {
			std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
			return -1;
		}
		std::string line, key, sep, val;
		while (std::getline(file_in, line))
		{
			std::stringstream ss(line);
			ss >> key; ss >> sep; ss >> val;
			inputparams_f[key] = val;
		}
		file_in.close();
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
				   csys = "PbPb";    if (inputparams_f.count("collsys")    > 0) csys    =      inputparams_f["collsys"];
									 if (  inputparams.count("collsys")    > 0) csys    =        inputparams["collsys"];
	                snn = "5020GeV"; if (inputparams_f.count("sNN")        > 0) snn     =      inputparams_f["sNN"];
						             if (  inputparams.count("sNN")        > 0) snn     =        inputparams["sNN"];

	              pname = "Charm";   if (inputparams_f.count("pName")      > 0) pname   =      inputparams_f["pName"];
						             if (  inputparams.count("pName")      > 0) pname   =        inputparams["pName"];

	               cent = "30-40%";  if (inputparams_f.count("centrality") > 0) cent    =      inputparams_f["centrality"];
						             if (  inputparams.count("centrality") > 0) cent    =        inputparams["centrality"];

	                 xb = 0.6;       if (inputparams_f.count("xB") 		   > 0) xb 	    = stod(inputparams_f["xB"]);
						             if (  inputparams.count("xB") 		   > 0) xb 	    = stod(  inputparams["xB"]);

	             eventn = 1000;      if (inputparams_f.count("eventN") 	   > 0) eventn  = stoi(inputparams_f["eventN"]);
						             if (  inputparams.count("eventN") 	   > 0) eventn  = stoi(  inputparams["eventN"]);

	std::string bcppstr = "20%";     if (inputparams_f.count("BCPP") 	   > 0) bcppstr =      inputparams_f["BCPP"];
						             if (  inputparams.count("BCPP") 	   > 0) bcppstr =        inputparams["BCPP"];

	            phiptsn = 25;        if (inputparams_f.count("phiGridN")   > 0) phiptsn = stoi(inputparams_f["phiGridN"]);
						             if (  inputparams.count("phiGridN")   > 0) phiptsn = stoi(  inputparams["phiGridN"]);

	              tstep = 0.1;       if (inputparams_f.count("TIMESTEP")   > 0) tstep   = stod(inputparams_f["TIMESTEP"]);
						             if (  inputparams.count("TIMESTEP")   > 0) tstep   = stod(  inputparams["TIMESTEP"]);

	              tcrit = 0.155;     if (inputparams_f.count("TCRIT") 	   > 0) tcrit   = stod(inputparams_f["TCRIT"]);
						             if (  inputparams.count("TCRIT") 	   > 0) tcrit   = stod(  inputparams["TCRIT"]);

	            bcpseed = 0;         if (inputparams_f.count("BCPSEED")    > 0) bcpseed = stoi(inputparams_f["BCPSEED"]);
						             if (  inputparams.count("BCPSEED")    > 0) bcpseed = stoi(  inputparams["BCPSEED"]);

	//checking if provided value of sNN is an option:
	if ((snn != "5440GeV") && (snn != "5020GeV") && (snn != "2760GeV") && (snn != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		return -2;
	}

	//setting BCPP double value:
	bcppstr.replace(bcppstr.find("%"), 1, ""); bcpp = stod(bcppstr)/100.0;

	return 1;
}

int GetInputs(int argc, char const *argv[], std::string &snn, std::string &pname, double &xb, size_t &LdndxMaxPts, size_t &LCollMaxPts, double &tcrit)
{

	std::vector<std::string> inputs; for (int i=2; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: -sNN=5020GeV -pName=Charm -xB=0.6 -LdndxMaxPoints=500000 -LCollMaxPoints=10000 --TCRIT=0.155" << std::endl;
		return 0;
	}

	std::map<std::string, std::string> inputparams;
	for (const auto &in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

	//checking if configuration file is provided:
	std::map<std::string, std::string> inputparams_f;
	if (inputparams.count("c") > 0) {
		std::ifstream file_in(inputparams["c"]);
		if (!file_in.is_open()) {
			std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
			return -1;
		}
		std::string line, key, sep, val;
		while (std::getline(file_in, line))
		{
			std::stringstream ss(line);
			ss >> key; ss >>sep; ss >> val;
			inputparams_f[key] = val;
		}
		file_in.close();
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
	        snn = "5020GeV"; if (inputparams_f.count("sNN")            > 0) snn         =      inputparams_f["sNN"];
						     if (inputparams.count("sNN")              > 0) snn         =        inputparams["sNN"];

	      pname = "Charm";   if (inputparams_f.count("pName")          > 0) pname       =      inputparams_f["pName"];
						     if (inputparams.count("pName")            > 0) pname       =        inputparams["pName"];

	         xb = 0.6;       if (inputparams_f.count("xB") 		       > 0) xb 	        = stod(inputparams_f["xB"]);
						     if (inputparams.count("xB") 		       > 0) xb 	        = stod(  inputparams["xB"]);

	LdndxMaxPts = 500000;    if (inputparams_f.count("LdndxMaxPoints") > 0) LdndxMaxPts = stoi(inputparams_f["LdndxMaxPoints"]);
						     if (  inputparams.count("LdndxMaxPoints") > 0) LdndxMaxPts = stoi(  inputparams["LdndxMaxPoints"]);

	LCollMaxPts = 10000;     if (inputparams_f.count("LCollMaxPoints") > 0) LCollMaxPts = stoi(inputparams_f["LCollMaxPoints"]);
						     if (inputparams.count("LCollMaxPoints")   > 0) LCollMaxPts = stoi(  inputparams["LCollMaxPoints"]);
	
		  tcrit = 0.155;     if (inputparams_f.count("TCRIT") 	   	   > 0) tcrit   	= stod(inputparams_f["TCRIT"]);
						    if (  inputparams.count("TCRIT") 	   	   > 0) tcrit    	= stod(  inputparams["TCRIT"]);

	//checking if provided value of sNN is an option:
	if ((snn != "5440GeV") && (snn != "5020GeV") && (snn != "2760GeV") && (snn != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		return -2;
	}

	return 1;
}