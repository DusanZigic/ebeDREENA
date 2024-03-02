/*//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

															ebeDREENA

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
#include "inputparser.hpp"
#include "mainheader.hpp"

#include <iostream>
#include <string>

int main(int argc, char const *argv[])
{
    if (argc == 1) {
		std::cout << "use -h for help" << std::endl;
		return 0;
	}

    if (argv[1] ==  std::string("-h")) {

		std::cout << "posible calculations: AverageEL, LTables" << std::endl;
		std::cout << "use: 'calculation' -h for extra help" << std::endl;
		return 0;
	}
    else if (argv[1] ==  std::string("AverageEL")) {
        // if (GetInputs(argc, argv, collsys, sNN, pName, centrality, xB, eventN, BCPP, phiGridN, TIMESTEP, TCRIT, BCPSEED) == 0) return -1;
        // AverageEL();
    }
    else if (argv[1] == std::string("LTables")) {

		// if (GetInputs(argc, argv, LT_sNN, LT_pName, LT_xB, LdndxMaxPoints, LCollMaxPoints) == 0) return -2;
		// GenerateLTables();
	}
	else {
		std::cerr << "posible calculations: AverageEL, LTables" << std::endl;
		return -3;
	}
}