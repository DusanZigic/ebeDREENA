#include "energyloss.hpp"

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

#include "utils.hpp"
#include "linearinterpolator.hpp"

void EnergyLoss::FdAHaltonSeqInit(std::size_t FdAMaxPts)
{
	m_FdAMaxPoints2 = FdAMaxPts;
	m_FdAMaxPoints3 = FdAMaxPts - 25;
	m_FdAMaxPoints4 = FdAMaxPts - 50;
	m_FdAMaxPoints5 = FdAMaxPts - 75; // NOTE (dusan): each consequent integral requires lesser precision
	
	m_FdAHS2.reserve(m_FdAMaxPoints2);
	m_FdAHS3.reserve(m_FdAMaxPoints3);
	m_FdAHS4.reserve(m_FdAMaxPoints4);
	m_FdAHS5.reserve(m_FdAMaxPoints5);
	
	for (std::size_t i = 0; i < FdAMaxPts; ++i) {
		m_FdAHS2.push_back(utils::haltonSequence((i + 1) * 409, 2));
		m_FdAHS3.push_back(utils::haltonSequence((i + 1) * 409, 3));
		m_FdAHS4.push_back(utils::haltonSequence((i + 1) * 409, 5));
		m_FdAHS5.push_back(utils::haltonSequence((i + 1) * 409, 7));
	}

	auto generateSortedIndices = [this](std::vector<std::size_t>& idxVec, std::size_t N) {
		idxVec.resize(N);
		std::iota(idxVec.begin(), idxVec.end(), 0);
		std::sort(idxVec.begin(), idxVec.end(), [this](std::size_t i, std::size_t j) { return m_FdAHS2[i] < m_FdAHS2[j]; });
	};

	generateSortedIndices(m_FdAHSSortedIdx2, m_FdAMaxPoints2);
	generateSortedIndices(m_FdAHSSortedIdx3, m_FdAMaxPoints3);
	generateSortedIndices(m_FdAHSSortedIdx4, m_FdAMaxPoints4);
	generateSortedIndices(m_FdAHSSortedIdx5, m_FdAMaxPoints5);
}

double EnergyLoss::dAp410(double ph, const LinearInterpolator<double> &norm) const noexcept
{
	return (
		1.0 / std::exp(norm.interpolate(ph))
	);
}

double EnergyLoss::FdA411(double ph, double dp, double invExpNorm, double ph_over_p, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * invExpNorm = 1.0 / std::exp(norm.interpolate(p))
// * ph_over_p  = ph / p = ph / (ph + dp)
{
	return (
		invExpNorm *
		dndx.interpolate(ph + dp, 1.0 - ph_over_p)
	);
}

double EnergyLoss::FdA412(double ph, double dp, double mFactor, double invExpNorm, double ph_over_p, const LinearInterpolator<double> &dndx, std::size_t p_idx_dndx) const noexcept
// NOTE (dusan): factored out variables:
// * mFactor      = m_mgC / (p + e) =  m_mgC / (p + std::sqrt(m_MC_sq + p * p))
// * invExpNorm   = 1.0 / std::exp(norm.interpolate(p))
// * ph_over_p    = ph / p = ph / (ph + dp)
// * m_mgC_over_2 = m_mgC / 2.0
// * p_idx_dndx   = dndx.locateIndex(0, p)
{
	if (dp < 2.0 * m_mgC_over_2) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph_over_p - mFactor;
	const double yq = yh - yl;
	double y;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints2; ++i) {
		const std::size_t sorted_idx = m_FdAHSSortedIdx2[i];

		y = yl + m_FdAHS2[sorted_idx] * yq;
		
		sum += dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y);
	}

	return (sum * invExpNorm * m_dAPoissonFactors[2] * yq / static_cast<double>(m_FdAMaxPoints2));
}

double EnergyLoss::FdA413(double ph, double dp, double mFactor, double invExpNorm, double ph_over_p, const LinearInterpolator<double> &dndx, std::size_t p_idx_dndx) const noexcept
// NOTE (dusan): factored out variables:
// * mFactor      = m_mgC / (p + e) =  m_mgC / (p + std::sqrt(m_MC_sq + p * p))
// * invExpNorm   = 1.0 / std::exp(norm.interpolate(p))
// * ph_over_p    = ph / p = ph / (ph + dp)
// * m_mgC_over_2 = m_mgC / 2.0
// * p_idx_dndx   = dndx.locateIndex(0, p)
{
	if (dp < 3.0 * m_mgC_over_2) return 0.0;
	
	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph_over_p - 2.0 * mFactor;
	const double yq = yh - yl;
	double y;

	const double zl = mFactor;
	double zh, zq, z;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints3; ++i) {
		const std::size_t sorted_idx = m_FdAHSSortedIdx3[i];

		y = yl + m_FdAHS2[sorted_idx] * yq;
		
		zh = 1.0 - ph_over_p - y - mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[sorted_idx] * zq;
		
		sum += dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) * zq;
	}

	return (sum * invExpNorm * m_dAPoissonFactors[3] * yq / static_cast<double>(m_FdAMaxPoints3));
}

double EnergyLoss::FdA414(double ph, double dp, double mFactor, double invExpNorm, double ph_over_p, const LinearInterpolator<double> &dndx, std::size_t p_idx_dndx) const noexcept
// NOTE (dusan): factored out variables:
// * mFactor      = m_mgC / (p + e) =  m_mgC / (p + std::sqrt(m_MC_sq + p * p))
// * invExpNorm   = 1.0 / std::exp(norm.interpolate(p))
// * ph_over_p    = ph / p = ph / (ph + dp)
// * m_mgC_over_2 = m_mgC / 2.0
// * p_idx_dndx   = dndx.locateIndex(0, p)
{
	if (dp < 4.0 * m_mgC_over_2) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph_over_p - 3.0 * mFactor;
	const double yq = yh - yl;
	double y;

	const double zl = mFactor;
	double zh, zq, z;

	const double zzl = mFactor;
	double zzh, zzq, zz;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_FdAMaxPoints4; ++i) {
		const std::size_t sorted_idx = m_FdAHSSortedIdx4[i];

		y = yl + m_FdAHS2[sorted_idx] * yq;
		
		zh = 1.0 - ph_over_p - y - 2.0 * mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[sorted_idx]*zq;
		
		zzh = 1.0 - ph_over_p - y - z - mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[sorted_idx]*zzq;
		
		sum += dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   zq * zzq;
	}

	return (sum * invExpNorm * m_dAPoissonFactors[4] * yq / static_cast<double>(m_FdAMaxPoints4));
}

double EnergyLoss::FdA415(double ph, double dp, double mFactor, double invExpNorm, double ph_over_p, const LinearInterpolator<double> &dndx, std::size_t p_idx_dndx) const noexcept
// NOTE (dusan): factored out variables:
// * mFactor      = m_mgC / (p + e) =  m_mgC / (p + std::sqrt(m_MC_sq + p * p))
// * invExpNorm   = 1.0 / std::exp(norm.interpolate(p))
// * ph_over_p    = ph / p = ph / (ph + dp)
// * m_mgC_over_2 = m_mgC / 2.0
// * p_idx_dndx   = dndx.locateIndex(0, p)
{
	if (dp < 5.0 * m_mgC_over_2) return 0.0;

	const double p = ph + dp;

	const double yl = mFactor;
	const double yh = 1.0 - ph_over_p - 4.0 * mFactor;
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
		const std::size_t sorted_idx = m_FdAHSSortedIdx5[i];

		y = yl + m_FdAHS2[sorted_idx] * yq;
		
		zh = 1.0 - ph_over_p - y - 3.0 * mFactor;
		zq = zh - zl;
		z = zl + m_FdAHS3[sorted_idx] * zq;

		zzh = 1.0 - ph_over_p - y - z - 2.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_FdAHS4[sorted_idx] * zzq;

		zzzh = 1.0 - ph_over_p - y - z - zz - mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_FdAHS5[sorted_idx] * zzzq;

		sum += dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz - zzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzz) *
			   zq * zzq * zzzq;

	}

	return (sum * invExpNorm * m_dAPoissonFactors[5] * yq / static_cast<double>(m_FdAMaxPoints5));
}

double EnergyLoss::FdA(double ph, double dp, const LinearInterpolator<double> &currnorm, const LinearInterpolator<double> &currdndx) const noexcept
{
	// NOTE (dusan): factor out everything that doesn't depend on variable of integration
	const double      p          = ph + dp;
	const double      e          = std::sqrt(m_MC_sq + p * p);
	const double      mFactor    = m_mgC / (p + e);
	const double      invExpNorm = std::exp(-currnorm.interpolate(p));
	const double      ph_over_p  = ph / p;
	const std::size_t p_idx_dndx = currdndx.locateIndex(0, p);

	return (
		FdA411(ph, dp,          invExpNorm, ph_over_p, currdndx) +
		FdA412(ph, dp, mFactor, invExpNorm, ph_over_p, currdndx, p_idx_dndx) +
		FdA413(ph, dp, mFactor, invExpNorm, ph_over_p, currdndx, p_idx_dndx) +
		FdA414(ph, dp, mFactor, invExpNorm, ph_over_p, currdndx, p_idx_dndx) +
		FdA415(ph, dp, mFactor, invExpNorm, ph_over_p, currdndx, p_idx_dndx)
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

	m_dAHS1.reserve(dAMaxPts);
	m_dAHS2.reserve(dAMaxPts);
	m_dAHS3.reserve(dAMaxPts);
	m_dAHS4.reserve(dAMaxPts);
	m_dAHS5.reserve(dAMaxPts);
	m_dAHS6.reserve(dAMaxPts);
	m_dAHS7.reserve(dAMaxPts);

	for (std::size_t i = 0; i < dAMaxPts; ++i) {
		m_dAHS1.push_back(utils::haltonSequence((i + 1) * 409,  2));
		m_dAHS2.push_back(utils::haltonSequence((i + 1) * 409,  3));
		m_dAHS3.push_back(utils::haltonSequence((i + 1) * 409,  5));
		m_dAHS4.push_back(utils::haltonSequence((i + 1) * 409,  7));
		m_dAHS5.push_back(utils::haltonSequence((i + 1) * 409, 11));
		m_dAHS6.push_back(utils::haltonSequence((i + 1) * 409, 13));
		m_dAHS7.push_back(utils::haltonSequence((i + 1) * 409, 17));
	}

	auto generateSortedIndices = [this](std::vector<std::size_t>& idxVec, std::size_t N) {
		idxVec.resize(N);
		std::iota(idxVec.begin(), idxVec.end(), 0);
		std::sort(idxVec.begin(), idxVec.end(), [this](std::size_t i, std::size_t j) { return m_dAHS1[i] < m_dAHS1[j]; });
	};

	generateSortedIndices(m_dAHSSortedIdx1, m_dAMaxPoints1);
	generateSortedIndices(m_dAHSSortedIdx2, m_dAMaxPoints2);
	generateSortedIndices(m_dAHSSortedIdx3, m_dAMaxPoints3);
	generateSortedIndices(m_dAHSSortedIdx4, m_dAMaxPoints4);
	generateSortedIndices(m_dAHSSortedIdx5, m_dAMaxPoints5);
	generateSortedIndices(m_dAHSSortedIdx6, m_dAMaxPoints6);
	generateSortedIndices(m_dAHSSortedIdx7, m_dAMaxPoints7);
}

double EnergyLoss::dA410(double ph, const LinearInterpolator<double> &norm) const noexcept
{
	return (
		m_dsdpti2.interpolate(ph) /
		std::exp(norm.interpolate(ph))
	);
}

double EnergyLoss::dA411(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + m_mgC_over_2;
	const double pq = p2 - p1;
	double p;
	double ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints1; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx1[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;
		ph_over_p = ph / p;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p);
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints1));
}

double EnergyLoss::dA412(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 2.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints2; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx2[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;
		
		yl = mFactor;
		yh = 1.0 - ph_over_p - mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[2] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   yq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints2));
}

double EnergyLoss::dA413(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 3.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double zl, zh, zq, z;

	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints3; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx3[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;
		
		yl = mFactor;
		yh = 1.0 - ph_over_p - 2.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		zl = mFactor;
		zh = 1.0 - ph_over_p - y - mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[sorted_idx] * zq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);

		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[3] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   yq * zq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints3));
}

double EnergyLoss::dA414(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 4.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;
	
	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints4; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx4[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;
		
		yl = mFactor;
		yh = 1.0 - ph_over_p - 3.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		zl = mFactor;
		zh = 1.0 - ph_over_p - y - 2.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[sorted_idx] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph_over_p - y - z - mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[sorted_idx] * zzq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[4] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   yq * zq * zzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints4));
}

double EnergyLoss::dA415(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 5.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints5; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx5[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;
		
		yl = mFactor;
		yh = 1.0 - ph_over_p - 4.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		zl = mFactor;
		zh = 1.0 - ph_over_p - y - 3.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[sorted_idx] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph_over_p - y - z - 2.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[sorted_idx] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph_over_p - y - z - zz - mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[sorted_idx] * zzzq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[5] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz - zzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzz) *
			   yq * zq * zzq * zzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints5));
}

double EnergyLoss::dA416(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 6.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double zzzzl, zzzzh, zzzzq, zzzz;
	
	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints6; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx6[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;
		
		yl = mFactor;
		yh = 1.0 - ph_over_p - 5.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		zl = mFactor;
		zh = 1.0 - ph_over_p - y - 4.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[sorted_idx] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph/p - y - z - 3.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[sorted_idx] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph_over_p - y - z - zz - 2.0 * mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[sorted_idx] * zzzq;

		zzzzl = mFactor;
		zzzzh = 1.0 - ph_over_p - y - z - zz - zzz - mFactor;
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[sorted_idx] * zzzzq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[6] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz - zzz - zzzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzzz) *
			   yq * zq * zzq * zzzq * zzzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints6));
}

double EnergyLoss::dA417(double ph, double p2, const LinearInterpolator<double> &norm, const LinearInterpolator<double> &dndx) const noexcept
// NOTE (dusan): factored out variables:
// * p2           = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0))
// * m_mgC_over_2 = m_mgC / 2.0
{
	const double p1 = ph + 7.0 * m_mgC_over_2;
	const double pq = p2 - p1;
	double p;

	double yl, yh, yq, y;
	
	double zl, zh, zq, z;
	
	double zzl, zzh, zzq, zz;
	
	double zzzl, zzzh, zzzq, zzz;
	
	double zzzzl, zzzzh, zzzzq, zzzz;
	
	double zzzzzl, zzzzzh, zzzzzq, zzzzz;
	
	double e, mFactor, ph_over_p;

	double sum = 0.0;
	for (std::size_t i = 0; i < m_dAMaxPoints7; ++i) {
		const std::size_t sorted_idx = m_dAHSSortedIdx7[i];

		p = p1 + m_dAHS1[sorted_idx] * pq;

		e = std::sqrt(m_MC_sq + p * p);
		mFactor = m_mgC / (p + e);
		ph_over_p = ph / p;

		yl = mFactor;
		yh = 1.0 - ph_over_p - 6.0 * mFactor;
		yq = yh - yl;
		y = yl + m_dAHS2[sorted_idx] * yq;

		zl = mFactor;
		zh = 1.0 - ph_over_p - y - 5.0 * mFactor;
		zq = zh - zl;
		z = zl + m_dAHS3[sorted_idx] * zq;

		zzl = mFactor;
		zzh = 1.0 - ph_over_p - y - z - 4.0 * mFactor;
		zzq = zzh - zzl;
		zz = zzl + m_dAHS4[sorted_idx] * zzq;

		zzzl = mFactor;
		zzzh = 1.0 - ph_over_p - y - z - zz - 3.0 * mFactor;
		zzzq = zzzh - zzzl;
		zzz = zzzl + m_dAHS5[sorted_idx] * zzzq;

		zzzzl = mFactor;
		zzzzh = 1.0 - ph_over_p - y - z - zz - zzz - 2.0 * mFactor;
		zzzzq = zzzzh - zzzzl;
		zzzz = zzzzl + m_dAHS6[sorted_idx] * zzzzq;

		zzzzzl = mFactor;
		zzzzzh = 1.0 - ph/p - y - z - zz - zzz - zzzz - mFactor;
		zzzzzq = zzzzzh - zzzzzl;
		zzzzz = zzzzzl + m_dAHS7[sorted_idx] * zzzzzq;

		std::size_t p_idx_dndx = dndx.locateIndex(0, p);
		
		sum += m_dsdpti2.interpolate(p) / p * std::exp(-norm.interpolate(p)) *
			   m_dAPoissonFactors[7] *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, 1.0 - ph_over_p - y - z - zz - zzz - zzzz - zzzzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, y) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, z) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzzz) *
			   dndx.interpolate_with_cached_x1(p_idx_dndx, p, zzzzz) *
			   yq* zq* zzq* zzzq* zzzzq* zzzzzq;
	}

	return (sum * pq / static_cast<double>(m_dAMaxPoints7));
}

double EnergyLoss::dA41(double ph, LinearInterpolator<double> &currnorm, LinearInterpolator<double> &currdndx) const noexcept
{
	// NOTE (dusan): factor out upper integration limit since it's contant for all dA integrals unlike lower limit
	const double p2 = (((2.0 * ph) < (ph + 30.0)) ? (2.0 * ph) : (ph + 30.0));

	if (m_pName == "Gluon") { // NOTE (dusan): gluons need 7 dA integrals
		return (
			dA410(ph,     currnorm) +
			dA411(ph, p2, currnorm, currdndx) +
			dA412(ph, p2, currnorm, currdndx) +
			dA413(ph, p2, currnorm, currdndx) +
			dA414(ph, p2, currnorm, currdndx) +
			dA415(ph, p2, currnorm, currdndx) +
			dA416(ph, p2, currnorm, currdndx) +
			dA417(ph, p2, currnorm, currdndx)
		);
	}
	else { // NOTE (dusan): light quarks need 5 dA integrals
		return (
			dA410(ph,     currnorm) +
			dA411(ph, p2, currnorm, currdndx) +
			dA412(ph, p2, currnorm, currdndx) +
			dA413(ph, p2, currnorm, currdndx) +
			dA414(ph, p2, currnorm, currdndx) +
			dA415(ph, p2, currnorm, currdndx)
		);
	}
}