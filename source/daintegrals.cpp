#include "energyloss.hpp"

#include <vector>
#include <cmath>

#include "utils.hpp"
#include "linearinterpolator.hpp"

void EnergyLoss::FdAHaltonSeqInit(std::size_t FdAMaxPts)
{
	m_FdAMaxPoints2 = FdAMaxPts;
	m_FdAMaxPoints3 = FdAMaxPts - 25;
	m_FdAMaxPoints4 = FdAMaxPts - 50;
	m_FdAMaxPoints5 = FdAMaxPts - 75; // NOTE (dusan): each consequent integral requires lesser precision
	
	for (std::size_t i=0; i<FdAMaxPts; i++)
	{
		m_FdAHS2.push_back(utils::haltonSequence((i + 1) * 409, 2));
		m_FdAHS3.push_back(utils::haltonSequence((i + 1) * 409, 3));
		m_FdAHS4.push_back(utils::haltonSequence((i + 1) * 409, 5));
		m_FdAHS5.push_back(utils::haltonSequence((i + 1) * 409, 7));
	}
}

double EnergyLoss::dAp410(double ph, const LinearInterpolator<double> &norm) const noexcept
{
	return (
		1.0 / std::exp(norm.interpolate(ph))
	);
}

double EnergyLoss::FdA411(double ph, double dp, double invExpNorm, const LinearInterpolator<double> &dndx) const noexcept
{
	return (
		invExpNorm *
		dndx.interpolate(ph + dp, 1.0 - ph / (ph + dp))
	);
}

double EnergyLoss::FdA412(double ph, double dp, double mFactor, double invExpNorm, const LinearInterpolator<double> &dndx) const noexcept
{
	if (dp < 2.0 * m_mgC / 2.0) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph / p - mFactor;
	const double yq = yh - yl;
	double y;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints2; ++i) {
		y = yl + m_FdAHS2[i] * yq;
		
		sum += dndx.interpolate(p, 1.0 - ph / p - y) *
			   dndx.interpolate(p, y);
	}

	return (sum * invExpNorm * m_dAPoissonFactors[2] * yq / static_cast<double>(m_FdAMaxPoints2));
}

double EnergyLoss::FdA413(double ph, double dp, double mFactor, double invExpNorm, const LinearInterpolator<double> &dndx) const noexcept
{
	if (dp < 3.0 * m_mgC / 2.0) return 0.0;
	
	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph / p - 2.0 * mFactor;
	const double yq = yh - yl;
	double y;

	const double zl = mFactor;
	double zh, zq, z;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints3; ++i) {
		y = yl + m_FdAHS2[i] * yq;
		
		zh = 1.0 - ph / p - y - mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[i] * zq;
		
		sum += dndx.interpolate(p, 1.0 - ph / p - y - z) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) * zq;
	}

	return (sum * invExpNorm * m_dAPoissonFactors[3] * yq / static_cast<double>(m_FdAMaxPoints3));
}

double EnergyLoss::FdA414(double ph, double dp, double mFactor, double invExpNorm, const LinearInterpolator<double> &dndx) const noexcept
{
	if (dp < 4.0*m_mgC / 2.0) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph / p - 3.0 * mFactor;
	const double yq = yh - yl;
	double y;

	const double zl = mFactor;
	double zh, zq, z;

	const double zzl = mFactor;
	double zzh, zzq, zz;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints4; ++i) {
		y = yl + m_FdAHS2[i] * yq;
		
		zh = 1.0 - ph / p - y - 2.0 * mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[i]*zq;
		
		zzh = 1.0 - ph / p - y - z - mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[i]*zzq;
		
		sum += dndx.interpolate(p, 1.0 - ph / p - y - z - zz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   zq * zzq;
	}

	return (sum * invExpNorm * m_dAPoissonFactors[4] * yq / static_cast<double>(m_FdAMaxPoints4));
}

double EnergyLoss::FdA415(double ph, double dp, double mFactor, double invExpNorm, const LinearInterpolator<double> &dndx) const noexcept
{
	if (dp < 5.0*m_mgC / 2.0) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph / p - 4.0 * mFactor;
	const double yq = yh - yl;
	double y;
		
	const double zl = mFactor;
	double zh, zq, z;

	const double zzl = mFactor;
	double zzh, zzq, zz;

	const double zzzl = mFactor;
	double zzzh, zzzq, zzz;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints5; ++i) {
		y = yl + m_FdAHS2[i] * yq;
		
		zh = 1.0 - ph / p - y - 3.0 * mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[i] * zq;

		zzh = 1.0 - ph / p - y - z - 2.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[i] * zzq;

		zzzh = 1.0 - ph / p - y - z - zz - mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_FdAHS5[i] * zzzq;

		sum += dndx.interpolate(p, 1.0 - ph / p - y - z - zz - zzz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   dndx.interpolate(p, zzz) *
			   zq * zzq * zzzq;

	}

	return (sum * invExpNorm * m_dAPoissonFactors[5] * yq / static_cast<double>(m_FdAMaxPoints5));
}

double EnergyLoss::FdA(double ph, double dp, const LinearInterpolator<double> &currnorm, const LinearInterpolator<double> &currdndx) const noexcept
{
	const double p          = ph + dp;
	const double e          = std::sqrt(m_MC_sq + p * p);
	const double mFactor    = m_mgC / (p + e);
	const double invExpNorm = 1.0 / std::exp(currnorm.interpolate(p)); // NOTE (dusan): factor out everything that doesn't depend on variable of integration

	return (
		FdA411(ph, dp,          invExpNorm, currdndx) +
		FdA412(ph, dp, mFactor, invExpNorm, currdndx) +
		FdA413(ph, dp, mFactor, invExpNorm, currdndx) +
		FdA414(ph, dp, mFactor, invExpNorm, currdndx) +
		FdA415(ph, dp, mFactor, invExpNorm, currdndx)
	);
}

void EnergyLoss::dAHaltonSeqInit(std::size_t dAMaxPts)
{
	m_dAMaxPoints1 = dAMaxPts;
	m_dAMaxPoints2 = dAMaxPts - 100;
	m_dAMaxPoints3 = dAMaxPts - 200;
	m_dAMaxPoints4 = dAMaxPts - 300;
	m_dAMaxPoints5 = dAMaxPts - 400;
	m_dAMaxPoints6 = dAMaxPts - 500;
	m_dAMaxPoints7 = dAMaxPts - 600; // NOTE (dusan): each consequent integral requires lesser precision

	for (std::size_t i=0; i<dAMaxPts; i++)
	{
		m_dAHS1.push_back(utils::haltonSequence((i + 1) * 409,  2));
		m_dAHS2.push_back(utils::haltonSequence((i + 1) * 409,  3));
		m_dAHS3.push_back(utils::haltonSequence((i + 1) * 409,  5));
		m_dAHS4.push_back(utils::haltonSequence((i + 1) * 409,  7));
		m_dAHS5.push_back(utils::haltonSequence((i + 1) * 409, 11));
		m_dAHS6.push_back(utils::haltonSequence((i + 1) * 409, 13));
		m_dAHS7.push_back(utils::haltonSequence((i + 1) * 409, 17));
	}
}

double EnergyLoss::dA410(double ph, const LinearInterpolator<double> &norm) const noexcept
{
	return (
		m_dsdpti2.interpolate(ph) /
		std::exp(norm.interpolate(ph))
	);
}

double EnergyLoss::dA411(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints1; ++i) {
		p = p1 + m_dAHS1[i] * pq;
		
		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   dndx.interpolate(p, 1.0 - ph / p);
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints1));
}

double EnergyLoss::dA412(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 2.0 * m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints2; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		
		yl = mFactor;
		yh = 1.0 - ph / p - mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;
		
		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[2] *
			   dndx.interpolate(p, 1.0 - ph / p - y) *
			   dndx.interpolate(p, y) *
			   yq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints2));
}

double EnergyLoss::dA413(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 3.0 * m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double zl, zh, zq, z;

	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints3; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		
		yl = mFactor;
		yh = 1.0 - ph / p - 2.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;

		zl = mFactor;
		zh = 1.0 - ph / p - y - mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[i] * zq;

		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[3] *
			   dndx.interpolate(p, 1.0 - ph / p - y - z) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   yq * zq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints3));
}

double EnergyLoss::dA414(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 4.0 * m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints4; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		
		yl = mFactor;
		yh = 1.0 - ph / p - 3.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;

		zl = mFactor;
		zh = 1.0 - ph / p - y - 2.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[i] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph / p - y - z - mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i] * zzq;

		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[4] *
			   dndx.interpolate(p, 1.0 - ph / p - y - z - zz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   yq * zq * zzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints4));
}

double EnergyLoss::dA415(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 5.0 * m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints5; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		
		yl = mFactor;
		yh = 1.0 - ph / p - 4.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;

		zl = mFactor;
		zh = 1.0 - ph / p - y - 3.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[i] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph / p - y - z - 2.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph / p - y - z - zz - mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i] * zzzq;

		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[5] *
			   dndx.interpolate(p, 1.0 - ph / p - y - z - zz - zzz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   dndx.interpolate(p, zzz) *
			   yq * zq * zzq * zzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints5));
}

double EnergyLoss::dA416(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 6.0*m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double zzzzl, zzzzh, zzzzq, zzzz;
	
	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints6; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		
		yl = mFactor;
		yh = 1.0 - ph / p - 5.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;

		zl = mFactor;
		zh = 1.0 - ph / p - y - 4.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[i] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph/p - y - z - 3.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph / p - y - z - zz - 2.0 * mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i] * zzzq;

		zzzzl = mFactor;
		zzzzh = 1.0 - ph / p - y - z - zz - zzz - mFactor;
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[i] * zzzzq;

		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[6] *
			   dndx.interpolate(p, 1.0 - ph / p - y - z - zz - zzz - zzzz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   dndx.interpolate(p, zzz) *
			   dndx.interpolate(p, zzzz) *
			   yq * zq * zzq * zzzq * zzzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints6));
}

double EnergyLoss::dA417(double ph, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
{
	const double p1 = ph + 7.0*m_mgC / 2.0;
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double zzzzl, zzzzh, zzzzq, zzzz;
	
	double zzzzzl, zzzzzh, zzzzzq, zzzzz;
	
	double e, mFactor;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints7; ++i) {
		p = p1 + m_dAHS1[i] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);

		yl = mFactor;
		yh = 1.0 - ph / p - 6.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[i] * yq;

		zl = mFactor;
		zh = 1.0 - ph / p - y - 5.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[i] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph / p - y - z - 4.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[i] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph / p - y - z - zz - 3.0 * mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[i] * zzzq;

		zzzzl = mFactor;
		zzzzh = 1.0 - ph / p - y - z - zz - zzz - 2.0 * mFactor;
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[i] * zzzzq;

		zzzzzl = mFactor;
		zzzzzh = 1.0 - ph/p - y - z - zz - zzz - zzzz - mFactor;
		zzzzzq = zzzzzh - zzzzzl;
		zzzzz = zzzzzl + m_dAHS7[i] * zzzzzq;

		sum += m_dsdpti2.interpolate(p) / p / std::exp(norm.interpolate(p)) *
			   m_dAPoissonFactors[7] *
			   dndx.interpolate(p, 1.0 - ph / p - y - z - zz - zzz - zzzz - zzzzz) *
			   dndx.interpolate(p, y) *
			   dndx.interpolate(p, z) *
			   dndx.interpolate(p, zz) *
			   dndx.interpolate(p, zzz) *
			   dndx.interpolate(p, zzzz) *
			   dndx.interpolate(p, zzzzz) *
			   yq* zq* zzq* zzzq* zzzzq* zzzzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints7));
}

double EnergyLoss::dA41(double ph, LinearInterpolator<double> &currnorm, LinearInterpolator<double> &currdndx) const noexcept
{
	if (m_pName == "Gluon") { // NOTE (dusan): gluons need 7 dA integrals
		return (
			dA410(ph, currnorm) +
			dA411(ph, currnorm, currdndx) +
			dA412(ph, currnorm, currdndx) +
			dA413(ph, currnorm, currdndx) +
			dA414(ph, currnorm, currdndx) +
			dA415(ph, currnorm, currdndx) +
			dA416(ph, currnorm, currdndx) +
			dA417(ph, currnorm, currdndx)
		);
	}
	else { // NOTE (dusan): light quarks need 5 dA integrals
		return (
			dA410(ph, currnorm) +
			dA411(ph, currnorm, currdndx) +
			dA412(ph, currnorm, currdndx) +
			dA413(ph, currnorm, currdndx) +
			dA414(ph, currnorm, currdndx) +
			dA415(ph, currnorm, currdndx)
		);
	}
}