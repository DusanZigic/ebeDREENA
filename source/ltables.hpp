#ifndef LTABLES_HPP
#define LTABLES_HPP

#include "utils.hpp"
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

    double m_xB_2;                               // squared xB (precalculated for optimizations)
    const double m_lambda_2 = m_lambda*m_lambda; // squared lambda (precalculated for optimizations)
    double m_alpha_prefactor;                    // prefactor for alpha (precalculated for optimizations)
    
    gridPoints m_Grids; //grids

    utils::ParticleType m_particleType;

    std::vector<double> m_LdndxHSeq1, m_LdndxHSeq2, m_LdndxHSeq3;
    void LdndxHSeqInit();
    
    std::vector<std::vector<std::vector<std::vector<double>>>> m_LdndxTbl;
    std::vector<std::vector<std::vector<double>>> m_LNormTbl;
    double dElossDYN(double tau, double x, double k, double q, double varphi, double T, double mu2, double mg2, double M2, double e, double b, double alpha1) const noexcept;
    double Ldndx(double tau, double T, double x, double mu2, double mg2, double M2, double e, double alpha1) const noexcept;
    void RadLTables();

    std::vector<double> m_LCollHSeq1, m_LCollHSeq2, m_LCollHSeq3;
    void LCollHSeqInit();
    
    std::vector<std::vector<double>> m_LCollTbl;
    std::complex<double> deltaL2(double q, double w, double mu2, double mu4) const noexcept;
    std::complex<double> deltaT2(double q, double w, double mu2, double mu4) const noexcept;
    double ENumFinite(double p, double T) const noexcept;
    void CollLTables();

    int exportLTables() const;

};

#endif // LTABLES_HPP