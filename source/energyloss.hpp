#ifndef HEADERFILE_ELOSSHEADER
#define HEADERFILE_ELOSSHEADER

#include "config.hpp"
#include "grids.hpp"
#include "linearinterpolator.hpp"

#include <string>
#include <vector>
#include <cstddef>

class EnergyLoss {

public:
    EnergyLoss(const config::energyLossConfig &cfg);
    ~EnergyLoss();
    void runEnergyLoss();

private:
    std::string m_collsys;	     // collision system
    std::string m_sNN; 		     // collision energy
    std::string m_pName; 	     // particle name
    std::string m_centrality;    // centrality class
    double m_xB;			     // xB value
    double m_BCPP;			     // binary collision points percentage
    std::size_t m_eventN;		 // number of events 
    std::size_t m_phiGridN;		 // phi points number
    double m_TIMESTEP;           // jets' traversal step in fm
    double m_TCRIT;	             // temperature at which eloss stops
    int m_BCPSEED;			     // seed for generating initial position points
    
    double m_nf;			     // effective number of flavours
    
    double m_mgC;                // constant gluon mass used for dA integrals
    double m_MC;	             // constant particle mass used for dA integrals
    double m_TCollConst;         // constant temperature used for Gauss filter integration
    
    gridPoints m_Grids;				    	              //grid points
    LinearInterpolator<double> m_LNorm, m_Ldndx, m_LColl; //interpolated L tables
    LinearInterpolator<double> m_dsdpti2; 				  //initial pT distribution    
    
    double m_tau0;									  	         //thermalization time
    double m_tauMaxFM = 25.0;					                 //maximal value of tau for TProfile grids in fm
    double m_tempTau0, m_tempTauStep;                            //temperature evolution tau grid parameters
    std::size_t m_tempTauMax;
    double m_tempX0, m_tempXStep;                                //temperature evolution x grid parameters
    std::size_t m_tempXMax;
    double m_tempY0, m_tempYStep;                                //temperature evolution y grid parameters
    std::size_t m_tempYMax;
    std::vector<double> m_tempTauGrid, m_tempXGrid, m_tempYGrid; //temperature evolutions grids
    
    std::vector<double> m_phiGridPts; //phi points

    std::size_t m_FdAMaxPoints2, m_FdAMaxPoints3, m_FdAMaxPoints4, m_FdAMaxPoints5; //number of points for FdA integration
    std::vector<double> m_FdAHS2, m_FdAHS3, m_FdAHS4, m_FdAHS5;                     //vectors that store Halton sequences for FdA integrals
    
    std::size_t m_dAMaxPoints1, m_dAMaxPoints2, m_dAMaxPoints3, m_dAMaxPoints4, m_dAMaxPoints5, m_dAMaxPoints6, m_dAMaxPoints7; //number of points for dA integration
    std::vector<double> m_dAHS1, m_dAHS2, m_dAHS3, m_dAHS4, m_dAHS5, m_dAHS6, m_dAHS7; 								 	        //vectors that store Halton sequences for dA integrals

    int loaddsdpti2(const std::string &pname, LinearInterpolator<double> &dsdpti2int) const;
    int loadLdndx();
    int loadLNorm();
    int loadLColl();
    int generateTempGrid();
    int loadPhiPoints();
    int loadBinCollPoints(std::size_t event_id, std::vector<std::vector<double>> &bcpoints);
    int generateInitPosPoints(std::size_t event_id, std::vector<double> &xPoints, std::vector<double> &yPoints);
    int loadTProfile(std::size_t event_id, LinearInterpolator<double> &tempProfile);

    void RadCollEL(double X0, double Y0, double phi0, const LinearInterpolator<double> &TProfile, std::vector<double> &radiativeRAA1, std::vector<std::vector<double>> &radiativeRAA2, std::vector<double> &collisionalEL, double &pathLength, double &temp) const;
    void RadCollEL(double X0, double Y0, double phi0, const LinearInterpolator<double> &TProfile, std::vector<double> &radiativeRAA, std::vector<double> &collisionalEL, double &pathLenght, double &temp) const;

    double haltonSequence(int index, int base) const;
    void FdAHaltonSeqInit(std::size_t FdAMaxPts);
    double dAp410(double ph, const LinearInterpolator<double> &norm) const;
    double FdA411(double ph, double dp, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double FdA412(double ph, double dp, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double FdA413(double ph, double dp, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double FdA414(double ph, double dp, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double FdA415(double ph, double dp, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double FdA(double ph, double dp, const LinearInterpolator<double> &currnorm, const LinearInterpolator<double> &currdndx) const;
    void dAHaltonSeqInit(std::size_t dAMaxPts);
    double dA410(double ph, const LinearInterpolator<double> &norm) const;
    double dA411(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA412(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA413(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA414(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA415(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA416(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA417(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const;
    double dA41(double ph, LinearInterpolator<double> &currnorm, LinearInterpolator<double> &currdndx) const;

    void generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab) const;
    void calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp) const;
    int exportResults(const std::string &particleName, std::size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, std::size_t trajecNum, std::size_t elossNum) const;

    void runELossHeavyFlavour();
    void gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const;

    void runELossLightQuarks();
    void gaussFilterIntegrate(const LinearInterpolator<double> &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const;

    void runELossLightFlavour();
    void gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA) const;
};

#endif