#ifndef LTABLES_HPP
#define LTABLES_HPP

#include "config.hpp"
#include "grids.hpp"

#include <string>
#include <vector>
#include <complex>
#include <cstddef>

class LTables {

public:
    LTables(const config::lTablesConfig &cfg);
    ~LTables();
    void runLTables();

private:
    std::string m_sNN;             // collision energy
    std::string m_pName;           // particle name
    double m_xB;                   // xB value
    std::size_t m_LdndxMaxPoints;  // number of qmc integration points
    std::size_t m_LCollMaxPoints;
    double m_TCRIT;                // critical temperature

    double m_nf;                   // effective number of flavours
    const double m_Ng = 3.0;	   // effective number of gluons
    const double m_lambda = 0.2;   // QCD scale
    const double m_kmaxColl = 5.0; // kMaxColl value
    double m_CR;		           // Casimir (3 for gluons, 4/3 for quakrs)
    
    gridPoints m_Grids; //grids

    double productLog(double x) const;
    double unitStep(double x) const;
    long double unitStep(long double x) const;

    std::vector<double> m_LdndxHSeq1, m_LdndxHSeq2, m_LdndxHSeq3;
    double haltonSequence(int index, int base) const;
    void LdndxHSeqInit();
    
    std::vector<std::vector<std::vector<std::vector<double>>>> m_LdndxTbl;
    std::vector<std::vector<std::vector<double>>> m_LNormTbl;
    double dElossDYN(double tau, double p, double x, double k, double q, double varphi, double T) const;
    double Ldndx(double tau, double p, double T, double x) const;
    void RadLTables();

    std::vector<double> m_LCollHSeq1, m_LCollHSeq2, m_LCollHSeq3;
    void LCollHSeqInit();
    
    std::vector<std::vector<double>> m_LCollTbl;
    std::complex<double> deltaL2(double q, double w, double T) const;
    std::complex<double> deltaT2(double q, double w, double T) const;
    double ENumFinite(double p, double T) const;
    void CollLTables();

    int exportLTables() const;

};

#endif // LTABLES_HPP