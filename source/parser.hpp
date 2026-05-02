#ifndef PARSER_HPP
#define PARSER_HPP

#include "config.hpp"

namespace parser {
    config::energyLossConfig parseEnergyLossArgs(int argc, const char* argv[]);
    config::lTablesConfig parseLTableArgs(int args, const char* argv[]);
} // namespace parser

#endif // PARSER_HPP