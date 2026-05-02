#include "parser.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>

namespace {
    std::string cleanKey(std::string key) {
        key.erase(std::remove(key.begin(), key.end(), '-'), key.end());
        const char* whitespace = " \t\n\r\f\v";
        key.erase(0, key.find_first_not_of(whitespace));
        key.erase(key.find_last_not_of(whitespace) + 1);
        return key;
    }

    std::string cleanVal(std::string val) {
        const char* whitespace = " \t\n\r\f\v";
        val.erase(0, val.find_first_not_of(whitespace));
        val.erase(val.find_last_not_of(whitespace) + 1);
        return val;
    }

    void parseCMDToMap(int argc, const char* argv[], const std::set<std::string>& validKeys, std::map<std::string, std::string>& cmdMap) {
        for (int i = 2; i < argc; ++i) {
            std::string arg = argv[i];
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                std::string key = cleanKey(arg.substr(0, pos));
                if (validKeys.find(key) == validKeys.end()) {
                    std::cerr << "FATAL: Unknown argument '" << key << "'" << std::endl;
                    std::exit(1);
                }
                std::string val = cleanVal(arg.substr(pos + 1));
                cmdMap[key] = val;
            }
        }
    }

    void loadFileToMap(const std::string& path, std::map<std::string, std::string>& target) {
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "FATAL: Unable to open config file: " << path << std::endl;
            std::exit(1);
        }
        std::string line, key, sep, val;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t sepPos = line.find('=');
            if (sepPos != std::string::npos) {
                std::string key = cleanKey(line.substr(0, sepPos));
                std::string val = cleanVal(line.substr(sepPos + 1));
                target[key] = val;
            }
        }
    }

    unsigned int to_uint(const std::string& s) {
        long val = std::stol(s);
        if (val < 0) {
            throw std::out_of_range("Negative value for unsigned int parameter");
        }
        return static_cast<unsigned int>(val);
    }

    std::size_t to_size_t(const std::string& s) {
        long long val = std::stoull(s);
        if (val < 0) {
            throw std::out_of_range("Negative value provided for a size_t parameter");
        }
        return static_cast<std::size_t>(val);
    }

    double to_double(const std::string& s) {
        return std::stod(s);
    }
    
    template <typename T, typename Func>
    void safeAssign(const std::string& key, T& target, 
                    const std::map<std::string, std::string>& fileMap,
                    const std::map<std::string, std::string>& cmdMap,
                    Func parser) {
        
        std::string source;
        bool found = false;
    
        if (fileMap.count(key)) { source = fileMap.at(key); found = true; }
        if (cmdMap.count(key))  { source = cmdMap.at(key); found = true; }
    
        if (found) {
            try {
                target = parser(source);
            } catch (const std::invalid_argument&) {
                std::cerr << "FATAL: Parameter '" << key << "' expects a numeric value, but got '" << source << "'" << std::endl;
                std::exit(1);
            } catch (const std::out_of_range&) {
                std::cerr << "FATAL: Parameter '" << key << "' value '" << source << "' is out of range for the data type." << std::endl;
                std::exit(1);
            }
        }
    }
}

namespace parser {
    config::lTablesConfig parseLTableArgs(int argc, const char* argv[]) {
        if (argc >= 3 && ((std::string(argv[2]) == "-h") || (std::string(argv[2]) == "--help"))) {
            config::lTablesConfig defaults;
            std::cout << "\n--- LTables Calculation Help ---\n"
                      << "Usage: " << argv[0] << " lTables [options]\n\n"
                      << "FORMATTING:\n"
                      << "  - Command line: Use '--key=value' (e.g., --pName=Charm)\n"
                      << "  - Config file:  Use 'key = value' (one per line)\n\n"
                      << "AVAILABLE PARAMETERS [Default Values]:\n"
                      << "  --sNN            [" << defaults.sNN << "]\n"
                      << "  --pName          [" << defaults.pName << "]\n"
                      << "  --xB             [" << defaults.xB << "]\n"
                      << "  --LdndxMaxPoints [" << defaults.LdndxMaxPoints << "]\n"
                      << "  --LCollMaxPoints [" << defaults.LCollMaxPoints << "]\n"
                      << "  --TCRIT          [" << defaults.TCRIT << "]\n"
                      << "\nOPTIONS:\n"
                      << "  --config=[path]  Load parameters from a file (overwritten by CLI)\n"
                      << "  -h, --help       Show this help message and exit\n" << std::endl;
            std::exit(0);
        }

        config::lTablesConfig cfg;
        std::map<std::string, std::string> cmdMap;
        std::set<std::string> validKeys = {
            "sNN", "pName", "xB", "LdndxMaxPoints", "LCollMaxPoints", "TCRIT", "config"
        };
        
        parseCMDToMap(argc, argv, validKeys, cmdMap);

        std::map<std::string, std::string> fileMap;
        if (cmdMap.count("config")) {
            loadFileToMap(cmdMap["config"], fileMap);
        }

        safeAssign("sNN",            cfg.sNN,              fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("pName",          cfg.pName,            fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("xB",             cfg.xB,               fileMap, cmdMap, to_double);
        safeAssign("LdndxMaxPoints", cfg.LdndxMaxPoints,   fileMap, cmdMap, to_size_t);
        safeAssign("LCollMaxPoints", cfg.LCollMaxPoints,   fileMap, cmdMap, to_size_t);
        safeAssign("TCRIT",          cfg.TCRIT,            fileMap, cmdMap, to_double);

        return cfg;
    }

    config::energyLossConfig parseEnergyLossArgs(int argc, const char* argv[]) {        
        if (argc >= 3 && ((std::string(argv[2]) == "-h") || (std::string(argv[2]) == "--help"))) {
            config::energyLossConfig defaults;
            std::cout << "\n--- Energy Loss Simulation Help ---\n"
                      << "Usage: " << argv[0] << " eLoss [options]\n\n"
                      << "FORMATTING:\n"
                      << "  - Command line: Use '--key=value' (e.g., --eventN=5000)\n"
                      << "  - Config file:  Use 'key = value' (one per line)\n\n"
                      << "AVAILABLE PARAMETERS [Default Values]:\n"
                      << "  --collsys    [" << defaults.collsys << "]\n"
                      << "  --sNN        [" << defaults.sNN << "]\n"
                      << "  --pName      [" << defaults.pName << "]\n"
                      << "  --centrality [" << defaults.centrality << "]\n"
                      << "  --xB         [" << defaults.xB << "]\n"
                      << "  --eventN     [" << defaults.eventN << "]\n"
                      << "  --BCPP       [" << defaults.BCPP <<"]\n"
                      << "  --phiGridN   [" << defaults.phiGridN << "]\n"
                      << "  --TIMESTEP   [" << defaults.TIMESTEP << "]\n"
                      << "  --TCRIT      [" << defaults.TCRIT << "]\n"
                      << "  --BCPSEED    [" << defaults.BCPSEED <<"]\n"
                      << "\nOPTIONS:\n"
                      << "  --config=[path]  Load parameters from a file (overwritten by CLI)\n"
                      << "  -h, --help       Show this help message and exit\n" << std::endl;
            std::exit(0);
        }
        
        config::energyLossConfig cfg;
        std::map<std::string, std::string> cmdMap;
        std::set<std::string> validKeys = {
            "collsys", "sNN", "pName", "centrality", "xB", 
            "eventN", "BCPP", "phiGridN", "TIMESTEP", "TCRIT", "BCPSEED", "config"
        };

        parseCMDToMap(argc, argv, validKeys, cmdMap);

        std::map<std::string, std::string> fileMap;
        if (cmdMap.count("config")) {
            loadFileToMap(cmdMap["config"], fileMap);
        }

        safeAssign("collsys",    cfg.collsys,    fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("sNN",        cfg.sNN,        fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("pName",      cfg.pName,      fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("centrality", cfg.centrality, fileMap, cmdMap, [](std::string v) { return v; });
        safeAssign("xB",         cfg.xB,         fileMap, cmdMap, to_double);
        safeAssign("eventN",     cfg.eventN,     fileMap, cmdMap, to_size_t);
        safeAssign("BCPP",       cfg.BCPP,       fileMap, cmdMap, to_double);
        safeAssign("phiGridN",   cfg.phiGridN,   fileMap, cmdMap, to_size_t);
        safeAssign("TIMESTEP",   cfg.TIMESTEP,   fileMap, cmdMap, to_double);
        safeAssign("TCRIT",      cfg.TCRIT,      fileMap, cmdMap, to_double);
        safeAssign("BCPSEED",    cfg.BCPSEED,    fileMap, cmdMap, to_uint);

        return cfg;
    }
} // namespace parser