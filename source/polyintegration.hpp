#ifndef HEADERFILE_POLYINTHEADER
#define HEADERFILE_POLYINTHEADER

#include <vector>

double linearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata);
double linearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);
double cubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata);
double cubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);

#endif