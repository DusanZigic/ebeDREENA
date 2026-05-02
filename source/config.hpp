#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <cstddef>

namespace config {
    struct lTablesConfig {
        std::string sNN = "5020GeV";         // collision energy
        std::string pName = "Charm";         // particle (capitalized)
        double xB = 0.6;                     // chromo -magnetic to -electric mass ratio
        std::size_t LdndxMaxPoints = 500000; // number of qmc integration points
        std::size_t LCollMaxPoints = 10000;
        double TCRIT = 0.155;                // temperature at which eloss stops
    };
    
    struct energyLossConfig {
        std::string collsys = "PbPb";      // collision system
        std::string sNN = "5020GeV";       // collision energy
        std::string pName = "Charm";       // particle (capitalized)
        std::string centrality = "30-40%"; // centrality class
        double xB = 0.6;                   // chromo -magnetic to -electric mass ratio
        std::size_t eventN = 1000;         // number of events
        double BCPP = 0.2;                 // ratio of binary collision points to be used as jets' initial positions
        std::size_t phiGridN = 25;         // number of jets' angles
        double TIMESTEP = 0.1;             // jets' traversal step in fm
        double TCRIT = 0.155;              // temperature at which eloss stops
        unsigned int BCPSEED = 0;          // seed for binary collision points shuffler (0 means no seed is set)
    };

} // namespace config

#endif // CONFIG_HPP