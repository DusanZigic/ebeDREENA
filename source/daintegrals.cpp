#include "elossheader.hpp"
#include "arsenal.hpp"
#include "linearinterpolation.hpp"

#include <vector>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//FdA INTEGRALS:

static size_t FdAMaxPoints2, FdAMaxPoints3, FdAMaxPoints4, FdAMaxPoints5; //number of points for dA integration
static std::vector<double> FdAHS2, FdAHS3, FdAHS4, FdAHS5;				  //vectors that store Halton sequences for dA integrals

//function that initializes Halton sequences for dA integrals calculations:
void FdAHaltonSeqInit(size_t FdAMaxPts)
{
	FdAMaxPoints2 = FdAMaxPts; 	  //setting values of dAMaxPoints
	FdAMaxPoints3 = FdAMaxPts-25;
	FdAMaxPoints4 = FdAMaxPts-50;
	FdAMaxPoints5 = FdAMaxPts-75;
	
	for (size_t i=0; i<FdAMaxPts; i++) //generating Halton sequences
	{
		FdAHS2.push_back(haltonSequence((i+1)*409, 2));
		FdAHS3.push_back(haltonSequence((i+1)*409, 3));
		FdAHS4.push_back(haltonSequence((i+1)*409, 5));
		FdAHS5.push_back(haltonSequence((i+1)*409, 7));
	}
}

double dAp410(double ph, interpolationF &normint) {
	return (1.0 / std::exp(normint.interpolation(ph)));
}

double FdA411(double ph, double dp, interpolationF &normint, interpolationF &dndxint) {
	return (1.0 / std::exp(normint.interpolation(ph + dp))*dndxint.interpolation(ph + dp, 1.0 - ph/(ph + dp)));
}

double FdA412(double ph, double dp, interpolationF &normint, interpolationF &dndxint) {
	if (dp < 2.0*mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double sum = 0.0;
	for (size_t i=0; i<FdAMaxPoints2; i++) {
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;
		sum += 1.0 / std::exp(normint.interpolation(p))*(1.0 / 2.0)*dndxint.interpolation(p, 1.0 - ph/p - y)*
			dndxint.interpolation(p, y)*(yh - yl);
	}

	return (sum/FdAMaxPoints2);
}

double FdA413(double ph, double dp, interpolationF &normint, interpolationF &dndxint) {
	if (dp < 3.0*mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double sum = 0.0;
	for (size_t i=0; i<FdAMaxPoints3; i++) {
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;
		sum += 1.0 / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum/FdAMaxPoints3);
}

double FdA414(double ph, double dp, interpolationF &normint, interpolationF &dndxint) {
	if (dp < 4.0*mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double sum = 0.0;
	for (size_t i=0; i<FdAMaxPoints4; i++) {
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + FdAHS4[i]*zzq;
		sum += 1.0 / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum/FdAMaxPoints4);
}

double FdA415(double ph, double dp, interpolationF &normint, interpolationF &dndxint) {
	if (dp < 5.0*mgC / 2.0) return 0.0;
	double p = ph + dp;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double sum = 0.0;
	for (size_t i=0; i<FdAMaxPoints5; i++) {
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 4.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + FdAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + FdAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + FdAHS4[i]*zzq;
		zzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - mgC/(p + std::sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + FdAHS5[i]*zzzq;
		sum += 1.0 / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*dndxint.interpolation(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);

	}

	return (sum/FdAMaxPoints5);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//dA INTEGRALS:

static size_t dAMaxPoints1, dAMaxPoints2, dAMaxPoints3, dAMaxPoints4, dAMaxPoints5, dAMaxPoints6, dAMaxPoints7; //number of points for dA integration
static std::vector<double> dAHS1, dAHS2, dAHS3, dAHS4, dAHS5, dAHS6, dAHS7; 								 	//vectors that store Halton sequences for dA integrals

//function that initializes Halton sequences for dA integrals calculations:
void dAHaltonSeqInit(size_t dAMaxPts)
{
	dAMaxPoints1 = dAMaxPts;	 //setting values of dAMaxPoints
	dAMaxPoints2 = dAMaxPts-100;
	dAMaxPoints3 = dAMaxPts-200;
	dAMaxPoints4 = dAMaxPts-300;
	dAMaxPoints5 = dAMaxPts-400;
	dAMaxPoints6 = dAMaxPts-500;
	dAMaxPoints7 = dAMaxPts-600;

	for (size_t i=0; i<dAMaxPts; i++) //generating Halton sequences
	{
		dAHS1.push_back(haltonSequence((i+1)*409, 2));
		dAHS2.push_back(haltonSequence((i+1)*409, 3));
		dAHS3.push_back(haltonSequence((i+1)*409, 5));
		dAHS4.push_back(haltonSequence((i+1)*409, 7));
		dAHS5.push_back(haltonSequence((i+1)*409, 11));
		dAHS6.push_back(haltonSequence((i+1)*409, 13));
		dAHS7.push_back(haltonSequence((i+1)*409, 17));
	}
}

double dA410(double ph, interpolationF &normint) {
	return (dsdpti2.interpolation(ph)/exp(normint.interpolation(ph)));
}

double dA411(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints1; i++) {
		p = p1 + dAHS1[i]*pq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*dndxint.interpolation(p, 1.0 - ph/p);
	}

	return (sum*pq/dAMaxPoints1);
}

double dA412(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + 2.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints2; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0)*dndxint.interpolation(p, 1.0 - ph / p - y)*
			dndxint.interpolation(p, y)*(yh - yl);
	}

	return (sum*pq / dAMaxPoints2);
}

double dA413(double ph, interpolationF &normint, interpolationF &dndxint)
{
	double p1 = ph + 3.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints3; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*(yh - yl)*(zh - zl);
	}

	return (sum*pq/dAMaxPoints3);
}

double dA414(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + 4.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints4; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*(yh - yl)*(zh - zl)*(zzh - zzl);
	}

	return (sum*pq/dAMaxPoints4);
}

double dA415(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + 5.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints5; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 4.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;
		zzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - mgC/(p + std::sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz - zzz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*dndxint.interpolation(p, zzz)*(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl);
	}

	return (sum*pq/dAMaxPoints5);
}

double dA416(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + 6.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double zzzzl, zzzzh, zzzzq, zzzz;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints6; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 5.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 4.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;
		zzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;
		zzzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + dAHS6[i]*zzzzq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz - zzz - zzzz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*dndxint.interpolation(p, zzz)*dndxint.interpolation(p, zzzz)*(yh - yl)*(zh - zl)*
			(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl);
	}

	return (sum*pq/dAMaxPoints6);
}

double dA417(double ph, interpolationF &normint, interpolationF &dndxint) {
	double p1 = ph + 7.0*mgC / 2.0;
	double p2 = (((2.0*ph) < (ph + 30.0)) ? (2.0*ph) : (ph + 30.0));
	double pq = p2 - p1;
	double p;
	double yl, yh, yq, y;
	double zl, zh, zq, z;
	double zzl, zzh, zzq, zz;
	double zzzl, zzzh, zzzq, zzz;
	double zzzzl, zzzzh, zzzzq, zzzz;
	double zzzzzl, zzzzzh, zzzzzq, zzzzz;
	double sum = 0.0;
	for (size_t i=0; i<dAMaxPoints7; i++) {
		p = p1 + dAHS1[i]*pq;
		yl = mgC/(p + std::sqrt(MC*MC + p*p));
		yh = 1.0 - ph/p - 6.0*mgC/(p + std::sqrt(MC*MC + p*p));
		yq = yh - yl;
		y = yl + dAHS2[i]*yq;
		zl = mgC/(p + std::sqrt(MC*MC + p*p));
		zh = 1.0 - ph/p - y - 5.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zq = zh - zl;
		z = zl + dAHS3[i]*zq;
		zzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzh = 1.0 - ph/p - y - z - 4.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzq = zzh - zzl;
		zz = zzl + dAHS4[i]*zzq;
		zzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzh = 1.0 - ph/p - y - z - zz - 3.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzzq = zzzh - zzzl;
		zzz = zzzl + dAHS5[i]*zzzq;
		zzzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzh = 1.0 - ph/p - y - z - zz - zzz - 2.0*mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + dAHS6[i]*zzzzq;
		zzzzzl = mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzzh = 1.0 - ph/p - y - z - zz - zzz - zzzz - mgC/(p + std::sqrt(MC*MC + p*p));
		zzzzzq = zzzzzh - zzzzzl;
		zzzzz = zzzzzl + dAHS7[i]*zzzzzq;
		sum += dsdpti2.interpolation(p) / p / std::exp(normint.interpolation(p))*(1.0 / 2.0 / 3.0 / 4.0 / 5.0 / 6.0 / 7.0)*dndxint.interpolation(p, 1.0 - ph/p - y - z - zz - zzz - zzzz - zzzzz)*
			dndxint.interpolation(p, y)*dndxint.interpolation(p, z)*dndxint.interpolation(p, zz)*dndxint.interpolation(p, zzz)*dndxint.interpolation(p, zzzz)*dndxint.interpolation(p, zzzzz)*
			(yh - yl)*(zh - zl)*(zzh - zzl)*(zzzh - zzzl)*(zzzzh - zzzzl)*(zzzzzh - zzzzzl);
	}

	return (sum*pq/dAMaxPoints7);
}