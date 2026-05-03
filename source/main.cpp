/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

															ebeDREENA

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
#include "config.hpp"
#include "parser.hpp"
#include "ltables.hpp"
#include "energyloss.hpp"

#include <iostream>

int main(int argc, const char* argv[])
{
	if (argc < 2) {
        parser::printGlobalUsage(argv[0]);
        return 1;
    }

	std::string mode = argv[1];

	if (mode == "lTables") {
		auto cfg = parser::parseLTableArgs(argc, argv);
		LTables lTables(cfg);
		// lTables.runLTables();
	} else if (mode == "eLoss") {
		auto cfg = parser::parseEnergyLossArgs(argc, argv);
		EnergyLoss energyLoss(cfg);
		// energyLoss.runEnergyLoss();
	} else {
		std::cerr << "Error: Unknown mode '" << mode << "'\n";
		parser::printGlobalUsage(argv[0]);
		return 1;
	}

	return 0;
}