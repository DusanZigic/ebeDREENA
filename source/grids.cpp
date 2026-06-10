#include "grids.hpp"
#include "utils.hpp"
#include "polyintegrator.hpp"
#include "linearinterpolator.hpp"

#include <vector>
#include <string>
#include <cmath>

gridPoints::gridPoints(const std::string &sNN, const std::string &particleName, double tcrit) {
	setGridPoints(sNN, particleName, tcrit);
}

gridPoints::gridPoints() {}

void gridPoints::setGridPointsBottom() {
	std::size_t taugridn = 21;
	std::vector<std::vector<double>> tauden{{0.0, 10.0}, {20.0, 10.0}};
	generateGrids(tauden, taugridn, m_tauPts);

	std::size_t pgridn = 25;
	double pgridmax = m_sNN == "200GeV" ? 100.0 : 200.0;
	std::vector<std::vector<double>> pden{{1.0, 8.0}, {20.0, 7.0}, {30.0, 3.0}, {60.0, 5.0}, {pgridmax, 1.0}};
	generateGrids(pden, pgridn, m_pPts);

	std::size_t Tgridn = 40;
	std::vector<std::vector<double>> Tden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(Tden, Tgridn, m_TPts);

	double mg = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_2;
	double M = 4.75;
	double MAXP = m_sNN == "200GeV" ? 100.0 : 200.0;
	std::size_t xgridn = 30;
	double xmin = mg/(MAXP + std::sqrt(MAXP*MAXP + M*M));
	for (std::size_t i=0; i<xgridn; i++)
		m_xPts.push_back(std::exp(std::log(xmin) - std::log(xmin)/static_cast<double>(xgridn-1)*static_cast<double>(i)));

	std::size_t Radgridn = 20;
	double Radgridmax = m_sNN == "200GeV" ? 70.0 : 170.0;
	std::vector<std::vector<double>> Radden{{2.0, 10.0}, {21.8, 10.0}, {44.5, 1.050001}, {Radgridmax, 1.0}};
	generateGrids(Radden, Radgridn, m_RadPts);

	double mgC = utils::debyeMass(m_nf, 3.0/2.0*m_TCRIT)*utils::constants::INV_SQRT_2;
	std::size_t Fdpgridn = 16;
	std::vector<std::vector<double>> Fdpden = {{5.0*mgC/2.0, 10.0}, {12.0, 5.0}, {30.0, 0.0}};
	generateGrids(Fdpden, Fdpgridn-4, m_FdpPts);
	m_FdpPts.insert(m_FdpPts.begin(), 4.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 3.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 2.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 1.0*mgC/2.0);

	std::size_t pCollgridn = 20;
	double pCollgridmax = m_sNN == "200GeV" ? 70.0 : 170.0;
	std::vector<std::vector<double>> pCollden{{1.0, 10.0}, {4.0, 10.0}, {9.0, 2.5}, {30.0, 0.6}, {60.0, 0.5}, {pCollgridmax, 0.3}};
	generateGrids(pCollden, pCollgridn, m_pCollPts);

	std::size_t TCollgridn = 40;
	std::vector<std::vector<double>> TCollden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(TCollden, TCollgridn, m_TCollPts);

	std::size_t fingridn = 30;
	double fingridmax = m_sNN == "200GeV" ? 50.0 : 150.0;
	std::vector<std::vector<double>> finden{{5.0, 10.0}, {50.0, 10.0}, {70.0, 5.0}, {fingridmax, 3.0}};
	generateGrids(finden, fingridn, m_finPts);
}

void gridPoints::setGridPointsCharm() {
	//tauPts:
	std::size_t taugridn = 21;
	std::vector<std::vector<double>> tauden{{0.0, 10.0}, {5.0, 10.0}, {10.0, 10.0}, {15.0, 10.0}, {20.0, 10.0}};
	generateGrids(tauden, taugridn, m_tauPts);
	
	//pPts:
	std::size_t pgridn = 25;
	double pgridmax = m_sNN == "200GeV" ? 100.0 : 200.0;
	std::vector<std::vector<double>> pden{{1.0, 8.0}, {20.0, 7.0}, {30.0, 3.0}, {60.0, 5.0}, {pgridmax, 1.0}};
	generateGrids(pden, pgridn, m_pPts);

	//TPts:
	std::size_t Tgridn = 40;
	std::vector<std::vector<double>> Tden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(Tden, Tgridn, m_TPts);

	//xPts:
	double mg = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_2;
	double M = 1.2;
	double MAXP = m_sNN == "200GeV" ? 100.0 : 200.0;
	std::size_t xgridn = 30;
	double xmin = mg/(MAXP + std::sqrt(MAXP*MAXP + M*M));
	for (std::size_t i=0; i<xgridn; i++)
		m_xPts.push_back(std::exp(std::log(xmin) - std::log(xmin)/static_cast<double>(xgridn-1)*static_cast<double>(i)));

	//RadPts:
	std::size_t Radgridn = 20;
	double Radgridmax = m_sNN == "200GeV" ? 70.0 : 170.0;
	std::vector<std::vector<double>> Radden{{2.0, 10.0}, {21.8, 10.0}, {44.5, 1.05}, {Radgridmax, 1.0}};
	generateGrids(Radden, Radgridn, m_RadPts);

	//FdpPts:
	double mgC = utils::debyeMass(m_nf, 3.0/2.0*m_TCRIT)*utils::constants::INV_SQRT_2;
	std::size_t Fdpgridn = 16;
	std::vector<std::vector<double>> Fdpden = {{5.0*mgC/2.0, 10.0}, {12.0, 5.0}, {30.0, 0.0}};
	generateGrids(Fdpden, Fdpgridn-4, m_FdpPts);
	m_FdpPts.insert(m_FdpPts.begin(), 4.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 3.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 2.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 1.0*mgC/2.0);

	//pCollPts:
	std::size_t pCollgridn = 20;
	double pCollgridmax = m_sNN == "200GeV" ? 70.0 : 170.0;
	std::vector<std::vector<double>> pCollden{{1.0, 10.0}, {4.0, 10.0}, {9.0, 2.5}, {30.0, 0.6}, {60.0, 0.5}, {pCollgridmax, 0.3}};
	generateGrids(pCollden, pCollgridn, m_pCollPts);

	//TCollPts:
	std::size_t TCollgridn = 40;
	std::vector<std::vector<double>> TCollden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(TCollden, TCollgridn, m_TCollPts);

	//finpts:
	std::size_t fingridn = 30;
	double fingridmax = m_sNN == "200GeV" ? 50.0 : 150.0;
	std::vector<std::vector<double>> finden{{5.0, 10.0}, {50.0, 10.0}, {70.0, 5.0}, {fingridmax, 3.0}};
	generateGrids(finden, fingridn, m_finPts);
}

void gridPoints::setGridPointsGluon() {
	//tauPts:
	std::size_t taugridn = 21;
	std::vector<std::vector<double>> tauden{{0.0, 10.0}, {20.0, 10.0}};
	generateGrids(tauden, taugridn, m_tauPts);

	//pPts:
	std::size_t pgridn = m_sNN == "200GeV" ? 35 : 50;
	double  pgridmax = m_sNN == "200GeV" ? 150.0 : 450.0;
	double pgridmaxw = m_sNN == "200GeV" ? 0.5 : 1.0;
	std::vector<std::vector<double>> pden{{1.0, 8.0}, {20.0, 7.0}, {40.0, 3.0}, {100.0, 5.0}, {pgridmax, pgridmaxw}};
	generateGrids(pden, pgridn, m_pPts);

	//TPts:
	std::size_t Tgridn = 40;
	std::vector<std::vector<double>> Tden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(Tden, Tgridn, m_TPts);

	//xPts:
	double mg = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_2;
	double  M = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_2;
	double MAXP = m_sNN == "200GeV" ? 150.0 : 450.0;
	std::size_t xgridn = 50;
	double xmin = mg/(MAXP + std::sqrt(MAXP*MAXP + M*M));
	for (std::size_t i=0; i<xgridn; i++)
		m_xPts.push_back(std::exp(std::log(xmin) - std::log(xmin)/static_cast<double>(xgridn-1)*static_cast<double>(i)));


	//RadPts:
	std::size_t Radgridn = m_sNN == "200GeV" ? 30 : 40;
	double Radgridmax = m_sNN == "200GeV" ? 120.0 : 420.0;
	std::vector<std::vector<double>> Radden{{2.0, 10.0}, {50.0, 10.0}, {70.0, 1.0}, {Radgridmax, 1.0}};
	generateGrids(Radden, Radgridn, m_RadPts);

	//FdpPts:
	double mgC = utils::debyeMass(m_nf, 3.0/2.0*m_TCRIT)*utils::constants::INV_SQRT_2;
	std::size_t Fdpgridn = 22;
	std::vector<std::vector<double>> Fdpden = {{5.0*mgC/2.0, 10.0}, {12.0, 5.0}, {30.0, 0.0}};
	generateGrids(Fdpden, Fdpgridn-4, m_FdpPts);
	m_FdpPts.insert(m_FdpPts.begin(), 4.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 3.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 2.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 1.0*mgC/2.0);

	//pCollPts:
	std::size_t pCollgridn = m_sNN == "200GeV" ? 30 : 40;
	double pCollgridmax = m_sNN == "200GeV" ? 120.0 : 420.0;
	std::vector<std::vector<double>> pCollden{{1.0, 10.0}, {4.0, 10.0}, {9.0, 2.5}, {30.0, 0.6}, {60.0, 0.5}, {pCollgridmax, 0.3}};
	generateGrids(pCollden, pCollgridn, m_pCollPts);

	//TCollPts:
	std::size_t TCollgridn = 40;
	std::vector<std::vector<double>> TCollden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(TCollden, TCollgridn, m_TCollPts);

	//finpts:
	std::size_t fingridn = m_sNN == "200GeV" ? 35 : 50;
	double fingridmax = m_sNN == "200GeV" ? 100.0 : 400.0;
	std::vector<std::vector<double>> finden{{5.0, 10.0}, {50.0, 10.0}, {70.0, 5.0}, {fingridmax, 3.0}};
	generateGrids(finden, fingridn, m_finPts);
}

void gridPoints::setGridPointsLQuarks() {
	std::size_t taugridn = 21;
	std::vector<std::vector<double>> tauden{{0.0, 10.0}, {20.0, 10.0}};
	generateGrids(tauden, taugridn, m_tauPts);

	std::size_t pgridn = m_sNN == "200GeV" ? 35 : 50;
	double  pgridmax = m_sNN == "200GeV" ? 150.0 : 450.0;
	double pgridmaxw = m_sNN == "200GeV" ? 0.5 : 1.0;
	std::vector<std::vector<double>> pden{{1.0, 8.0}, {20.0, 7.0}, {40.0, 3.0}, {100.0, 5.0}, {pgridmax, pgridmaxw}};
	generateGrids(pden, pgridn, m_pPts);

	std::size_t Tgridn = 40;
	std::vector<std::vector<double>> Tden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(Tden, Tgridn, m_TPts);

	double mg = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_2;
	double M  = utils::debyeMass(m_nf, m_TPts[0])*utils::constants::INV_SQRT_6;
	double MAXP = m_sNN == "200GeV" ? 150.0 : 450.0;
	std::size_t xgridn = 50;
	double xmin = mg/(MAXP + std::sqrt(MAXP*MAXP + M*M));
	for (std::size_t i=0; i<xgridn; i++)
		m_xPts.push_back(std::exp(std::log(xmin) - std::log(xmin)/static_cast<double>(xgridn-1)*static_cast<double>(i)));

	std::size_t Radgridn = m_sNN == "200GeV" ? 30 : 40;
	double Radgridmax = m_sNN == "200GeV" ? 120.0 : 420.0;
	std::vector<std::vector<double>> Radden{{2.0, 10.0}, {50.0, 10.0}, {70.0, 1.0}, {Radgridmax, 1.0}};
	generateGrids(Radden, Radgridn, m_RadPts);

	double mgC = utils::debyeMass(m_nf, 3.0/2.0*m_TCRIT)*utils::constants::INV_SQRT_2;
	std::size_t Fdpgridn = 22;
	std::vector<std::vector<double>> Fdpden = {{5.0*mgC/2.0, 10.0}, {12.0, 5.0}, {30, 0.0}};
	generateGrids(Fdpden, Fdpgridn-4, m_FdpPts);
	m_FdpPts.insert(m_FdpPts.begin(), 4.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 3.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 2.0*mgC/2.0);
	m_FdpPts.insert(m_FdpPts.begin(), 1.0*mgC/2.0);

	std::size_t pCollgridn = m_sNN == "200GeV" ? 30 : 40;
	double pCollgridmax = m_sNN == "200GeV" ? 120.0 : 420.0;
	std::vector<std::vector<double>> pCollden{{1.0, 10.0}, {4.0, 10.0}, {9.0, 2.5}, {30.0, 0.6}, {60.0, 0.5}, {pCollgridmax, 0.3}};
	generateGrids(pCollden, pCollgridn, m_pCollPts);

	std::size_t TCollgridn = 40;
	std::vector<std::vector<double>> TCollden{{0.01, 10.0}, {2.0, 10.0}};
	generateGrids(TCollden, TCollgridn, m_TCollPts);

	std::size_t fingridn = m_sNN == "200GeV" ? 35 : 50;
	double fingridmax = m_sNN == "200GeV" ? 100.0 : 400.0;
	std::vector<std::vector<double>> finden{{5.0, 10.0}, {50.0, 10.0}, {70.0, 5.0}, {fingridmax, 3.0}};
	generateGrids(finden, fingridn, m_finPts);
}

void gridPoints::roundGridPoints() {
	for (std::size_t i=0; i<m_tauPts.size(); i++) m_tauPts[i] = std::round(m_tauPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_pPts.size();   i++)   m_pPts[i] = std::round(m_pPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_TPts.size();   i++)   m_TPts[i] = std::round(m_TPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_xPts.size();   i++)   m_xPts[i] = std::round(m_xPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_RadPts.size(); i++) m_RadPts[i] = std::round(m_RadPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_FdpPts.size(); i++) m_FdpPts[i] = std::round(m_FdpPts[i]*1e10)/1e10;

	for (std::size_t i=0; i<m_pCollPts.size(); i++) m_pCollPts[i] = std::round(m_pCollPts[i]*1e10)/1e10;
	for (std::size_t i=0; i<m_TCollPts.size(); i++) m_TCollPts[i] = std::round(m_TCollPts[i]*1e10)/1e10;

	for (std::size_t i=0; i<m_finPts.size(); i++) m_finPts[i] = std::round(m_finPts[i]*1e10)/1e10;
}

void gridPoints::setGridPoints(const std::string &sNN, const std::string &particleName, double tcrit)
{
	m_sNN = sNN;
	m_nf = m_sNN == "200GeV" ? 2.5 : 3.0;
	m_TCRIT = tcrit;

	if (particleName == "Bottom") {
		setGridPointsBottom();		
	} else if (particleName == "Charm") {
		setGridPointsCharm();		
	} else if (particleName == "Gluon") {
		setGridPointsGluon();	
	} else {
		setGridPointsLQuarks();
	}

	roundGridPoints();
}

gridPoints::~gridPoints() {}

const std::vector<double> & gridPoints::tauPts() const {
    return m_tauPts;
}
double gridPoints::tauPts(int i) const {
    if (i < 0)
        return m_tauPts.at(m_tauPts.size() + i);
    return m_tauPts.at(i);
}
std::size_t gridPoints::tauPtsLength() const {
    return m_tauPts.size();
}

const std::vector<double> & gridPoints::pPts() const {
    return m_pPts;
}
double gridPoints::pPts(int i) const {
    if (i < 0) return m_pPts.at(m_pPts.size() + i);
    return m_pPts.at(i);
}
std::size_t gridPoints::pPtsLength() const {
    return m_pPts.size();
}

const std::vector<double> & gridPoints::TPts() const {
    return m_TPts;
}
double gridPoints::TPts(int i) const {
    if (i < 0)
        return m_TPts.at(m_TPts.size() + i);
    return m_TPts.at(i);
}
std::size_t gridPoints::TPtsLength() const {
    return m_TPts.size();
}

const std::vector<double> & gridPoints::xPts() const {
    return m_xPts;
}
double gridPoints::xPts(int i) const {
    if (i < 0)
        return m_xPts.at(m_xPts.size() + i);
    return m_xPts.at(i);
}
std::size_t gridPoints::xPtsLength() const {
    return m_xPts.size();
}

const std::vector<double> & gridPoints::RadPts() const {
    return m_RadPts;
}
double gridPoints::RadPts(int i) const {
    if (i < 0)
        return m_RadPts.at(m_RadPts.size() + i);
    return m_RadPts.at(i);
}
std::size_t gridPoints::RadPtsLength() const {
    return m_RadPts.size();
}

const std::vector<double> & gridPoints::FdpPts() const {
    return m_FdpPts;
}
double gridPoints::FdpPts(int i) const {
    if (i < 0)
        return m_FdpPts.at(m_FdpPts.size() + i);
    return m_FdpPts.at(i);
}
std::size_t gridPoints::FdpPtsLength() const {
    return m_FdpPts.size();
}

const std::vector<double> & gridPoints::pCollPts() const {
    return m_pCollPts;
}
double gridPoints::pCollPts(int i) const {
    if (i < 0)
        return m_pCollPts.at(m_pCollPts.size() + i);
    return m_pCollPts.at(i);
}
std::size_t gridPoints::pCollPtsLength() const {
    return m_pCollPts.size();
}

const std::vector<double> & gridPoints::TCollPts() const {
    return m_TCollPts;
}
double gridPoints::TCollPts(int i) const {
    if (i < 0)
        return m_TCollPts.at(m_TCollPts.size() + i);
    return m_TCollPts.at(i);
}
std::size_t gridPoints::TCollPtsLength() const {
    return m_TCollPts.size();
}

const std::vector<double> & gridPoints::finPts() const {
    return m_finPts;
}
double gridPoints::finPts(int i) const {
    if (i < 0)
        return m_finPts.at(m_finPts.size() + i);
    return m_finPts.at(i);
}
std::size_t gridPoints::finPtsLength() const {
    return m_finPts.size();
}

void gridPoints::generateGrids(const std::vector<std::vector<double>> &density, std::size_t numpts, std::vector<double> &gridpoints)
{
	std::vector<double> densityX, densityF;
	for (std::size_t i=0; i<density.size(); i++) {densityX.push_back(density[i][0]); densityF.push_back(density[i][1]);}

	std::vector<double> inttabX, inttabF;

	double xxx = densityX.front();
	inttabX.push_back(0.0);
	inttabF.push_back(xxx);

	for (std::size_t i=1; i<19; i++)
	{
		xxx = densityX.front() + (densityX.back()-densityX.front())/static_cast<double>(19)*static_cast<double>(i);
		inttabX.push_back(poly::linearIntegrate(densityX, densityF, densityX[0], xxx));
		inttabF.push_back(xxx);
	}

	xxx = densityX.back();
	inttabX.push_back(poly::linearIntegrate(densityX, densityF, densityX[0], xxx));
	inttabF.push_back(xxx);

	LinearInterpolator<double> inttabInt(inttabX, inttabF);

	gridpoints.resize(0);

	gridpoints.push_back(densityX.front());

	for (std::size_t i=1; i<numpts-1; i++)
	{
		double a = inttabX.front() + (inttabX.back()-inttabX.front())*static_cast<double>(i)/static_cast<double>(numpts-1);
		gridpoints.push_back(inttabInt.interpolate(a));
	}

	gridpoints.push_back(densityX.back());
}