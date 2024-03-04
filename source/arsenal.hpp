#ifndef HEADERFILE_ARSENALHEADER
#define HEADERFILE_ARSENALHEADER

#include "linearinterpolation.hpp"

double productLog(double x);
double unitStep(double x);
long double unitStep(long double x);

double haltonSequence(int index, int base);

int generateTempGrid();

void generatePhiGridPts();
int generateInitPosPoints(size_t event_id, std::vector<double> &xPoints, std::vector<double> &yPoints);

void gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2);                                             //function that performs Gauss filter integration - modefied pT integration algorithm
void gaussFilterIntegrate(const interpolationF &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2); //function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
void gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA);                                                                                                                                           //function that performs Gauss filter integration - default algorithm

void calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp); //function that calculates averaged path-length and temperatures

#endif