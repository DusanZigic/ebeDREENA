#ifndef HEADERFILE_GRIDPOINTS
#define HEADERFILE_GRIDPOINTS

#include <vector>
#include <string>
#include <cstddef>

class gridPoints {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//public functions:
public:

	//CONSTRUCTORS:
	gridPoints();
	gridPoints(const std::string &sNN, const std::string &particleName, double tcrit);
	void setGridPoints(const std::string &sNN, const std::string &particleName, double tcrit);

	//DESTRUCTOR:
	~gridPoints();

	//GRID FUNCTIONS:
	const std::vector<double> & tauPts() const;
	double tauPts(int i) const;
	std::size_t tauPtsLength() const;

	const std::vector<double> & pPts() const;
	double pPts(int i) const;
	std::size_t pPtsLength() const;

	const std::vector<double> & xPts() const;
	double xPts(int i) const;
	std::size_t xPtsLength() const;

	const std::vector<double> & TPts() const;
	double TPts(int i) const;
	std::size_t TPtsLength() const;

	const std::vector<double> & FdpPts() const;
	double FdpPts(int i) const;
	std::size_t FdpPtsLength() const;

	const std::vector<double> & RadPts() const;
	double RadPts(int i) const;
	std::size_t RadPtsLength() const;

	const std::vector<double> & pCollPts() const;
	double pCollPts(int i) const;
	std::size_t pCollPtsLength() const;

	const std::vector<double> & TCollPts() const;
	double TCollPts(int i) const;
	std::size_t TCollPtsLength() const;

	const std::vector<double> & finPts() const;
	double finPts(int i) const;
	std::size_t finPtsLength() const;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private variables and functions:
private:

	std::string m_sNN;
	double m_nf;
	double m_TCRIT;
	std::vector<double> m_tauPts, m_pPts, m_TPts, m_xPts, m_RadPts, m_FdpPts;
	std::vector<double> m_pCollPts, m_TCollPts, m_finPts;
	void generateGrids(const std::vector<std::vector<double>> &density, std::size_t numpts, std::vector<double> &gridpoints);
	void setGridPointsBottom();
	void setGridPointsCharm();
	void setGridPointsGluon();
	void setGridPointsLQuarks();
	void roundGridPoints();
};

#endif