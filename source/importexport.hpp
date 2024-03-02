#ifndef HEADERFILE_IMPORTEXPORTHEADER
#define HEADERFILE_IMPORTEXPORTHEADER

#include "linearinterpolation.hpp"

#include <string>
#include <vector>

int loaddsdpti2();
int loaddsdpti2(const std::string &pname, interpolationF &dsdpti2int);
int loadLdndx();
int loadLNorm();
int loadLColl();
int loadTempGridParams();
int loadTProfile(size_t event_id, interpolationF &tempProfile);
int loadBinCollPoints(size_t event_id, std::vector<std::vector<double>> &bcpoints);
int loadPhiPoints();

int exportResults(size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, size_t trajecNum, size_t elossNum);
int exportResults(const std::string &particleName, size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, size_t trajecNum, size_t elossNum);
int exportLTables(const std::vector<double> &ldndxTable, const std::vector<double> &lnormTable, const std::vector<double> &lcollTable);

#endif