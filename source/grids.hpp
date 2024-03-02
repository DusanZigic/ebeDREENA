#ifndef HEADERFILE_GRIDPOINTS
#define HEADERFILE_GRIDPOINTS

#include <vector>
#include <string>

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
	double tauPts(size_t i) const;
	size_t tauPtsLength() const;

	const std::vector<double> & pPts() const;
	double pPts(size_t i) const;
	size_t pPtsLength() const;

	const std::vector<double> & xPts() const;
	double xPts(size_t i) const;
	size_t xPtsLength() const;

	const std::vector<double> & TPts() const;
	double TPts(size_t i) const;
	size_t TPtsLength() const;

	const std::vector<double> & FdpPts() const;
	double FdpPts(size_t i) const;
	size_t FdpPtsLength() const;

	const std::vector<double> & RadPts() const;
	double RadPts(size_t i) const;
	size_t RadPtsLength() const;

	const std::vector<double> & pCollPts() const;
	double pCollPts(size_t i) const;
	size_t pCollPtsLength() const;

	const std::vector<double> & TCollPts() const;
	double TCollPts(size_t i) const;
	size_t TCollPtsLength() const;

	const std::vector<double> & finPts() const;
	double finPts(size_t i) const;
	size_t finPtsLength() const;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private variables and functions:
private:

	double m_nf     = 3.0;
	double m_lambda = 0.2;
	double m_TCRIT  = 0.155;
	double productLog(double x);
	double muF(double temp);
	std::vector<double> m_tauPts, m_pPts, m_TPts, m_xPts, m_RadPts, m_FdpPts;
	std::vector<double> m_pCollPts, m_TCollPts, m_finPts;
	double linearIntegrate(const std::vector<double> &dataX, const std::vector<double> &dataF, double xH) const;
	std::vector<double> generateGrids(const std::vector<std::vector<double>> &density, size_t numpts) const;
};

#endif