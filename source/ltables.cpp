#include "ltables.hpp"
#include "grids.hpp"
#include "polyintegrator.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>

LTables::LTables(const config::lTablesConfig &cfg) {
	m_sNN            = cfg.sNN;
    m_pName          = cfg.pName;
    m_xB             = cfg.xB;
    m_LdndxMaxPoints = cfg.LdndxMaxPoints;
    m_LCollMaxPoints = cfg.LCollMaxPoints;
    m_TCRIT          = cfg.TCRIT;

	m_nf = m_sNN   == "200GeV" ? 2.5 : 3.0;
	m_CR = m_pName == "Gluon"  ? 3.0 : 4.0/3.0;

	m_xB_sq = m_xB*m_xB;
	m_alpha_prefactor = 4.0*utils::constants::PI/(11.0 - 2.0*m_nf/3.0);

	if (m_pName == "Bottom") {
		m_particleType = utils::ParticleType::Bottom;
	} else if (m_pName == "Charm") {
		m_particleType = utils::ParticleType::Charm;
	} else if (m_pName == "Gluon") {
		m_particleType = utils::ParticleType::Gluon;
	} else {
		m_particleType = utils::ParticleType::LQuarks;
	}
}

LTables::~LTables() {}

void LTables::runLTables() {
	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

    RadLTables();

    CollLTables();

	if (exportLTables() != 1) return;
}

void LTables::LdndxHSeqInit() {
	m_LdndxHSeq1.resize(m_LdndxMaxPoints);
	m_LdndxHSeq2.resize(m_LdndxMaxPoints);
	m_LdndxHSeq3.resize(m_LdndxMaxPoints);
	for (std::size_t i = 0; i < m_LdndxMaxPoints; ++i) {
		m_LdndxHSeq1[i] = utils::haltonSequence((i+1)*409, 2);
		m_LdndxHSeq2[i] = utils::haltonSequence((i+1)*409, 3);
		m_LdndxHSeq3[i] = utils::haltonSequence((i+1)*409, 5);
	}
}

double LTables::dElossDYN(double tau, double x, double x_sq, double k, double q, double varphi, double T, double mu_sq, double mg_sq, double M_sq, double e, double b, double alpha1) const noexcept {	
	double k_sq = k*k;
	double q_sq = q*q;
	double b_sq = b*b;
	double k_q_cosvarphi = k*q*std::cos(varphi);
	double k_sq_q_sq_plus_2_k_q_cosvarphi = k_sq + q_sq + 2.0*k_q_cosvarphi;
	double k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq = k_sq_q_sq_plus_2_k_q_cosvarphi + b_sq;
	double k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq_squared = k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq*k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq;
	
	double alpha  = m_alpha_prefactor/std::log((k_sq + mg_sq + M_sq*x_sq)/x*utils::constants::INV_LAMBDA_QCD_SQ);
	
	double fn = 1.0;
	fn *= 1.0*utils::constants::INV_HBARC_GEVFM * m_CR*alpha*utils::constants::INV_PI * 3.0*alpha1*T * 2.0*k*q*utils::constants::INV_PI;
	fn *= (mu_sq - mu_sq*m_xB_sq)/(q_sq + mu_sq*m_xB_sq)/(q_sq + mu_sq);

	double psi = (k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq)/2.0/x/e*tau*utils::constants::INV_HBARC_GEVFM;

	fn *= (1 - std::cos(psi));
	fn *= 2.0/(k_sq + b_sq)/k_sq_q_sq_2_k_q_cosvarphi_plus_b_sq_squared;
	fn *= (-1.0*k_q_cosvarphi*k_sq_q_sq_plus_2_k_q_cosvarphi + b_sq*(k_q_cosvarphi + q_sq));

	return fn;
}

double LTables::Ldndx(double tau, double T, double x, double mu_sq, double mg_sq, double M_sq, double e, double alpha1) const noexcept {
	double x_sq = x*x;
	double b = std::sqrt(mg_sq + M_sq*x_sq);
	double kl = 0.00000001; 
	double kh = 2.0*x*(1 - x)*e;
	double kq = (kh - kl);
	double ql = 0.000001;
	double qh = std::sqrt(4.0*e*T);
	double qq = qh - ql;
	double phil = 0.0;
	double phih = utils::constants::PI;
	double phiq = (phih - phil);
	double sum = 0.0;
	double k, q, phi;

	#pragma omp parallel for reduction(+:sum) private(k,q,phi)
	for (std::size_t i = 0; i<m_LdndxMaxPoints; i++) {
		  k  =   kl + m_LdndxHSeq1[i]*kq;
		  q  =   ql + m_LdndxHSeq2[i]*qq;
		phi  = phil + m_LdndxHSeq3[i]*phiq;
		sum += 2.0*dElossDYN(tau, x, x_sq, k, q, phi, T, mu_sq, mg_sq, M_sq, e, b, alpha1) / x;
	}

	return (sum*kq*qq*phiq/static_cast<double>(m_LdndxMaxPoints));
}

void LTables::RadLTables() {
	LdndxHSeqInit();
	
	const std::vector<double>& tauPts = m_Grids.tauPts();
	const std::vector<double>& pPts = m_Grids.pPts();
	const std::vector<double>& TPts = m_Grids.TPts();
	const std::vector<double>& xPts = m_Grids.xPts();

	m_LdndxTbl.resize(
		tauPts.size(), std::vector<std::vector<std::vector<double>>>(
			pPts.size(), std::vector<std::vector<double>>(
				TPts.size(), std::vector<double>(
					xPts.size(), 0.0
				)
			)
		)
	);
	m_LNormTbl.resize(
		tauPts.size(), std::vector<std::vector<double>>(
			pPts.size(), std::vector<double>(
				TPts.size(), 0.0
			)
		)
	);
	
	double xIntegLimitLow, xIntegLimitHigh;

	for (std::size_t i_tau = 0; i_tau < tauPts.size(); ++i_tau) {
		for (std::size_t i_p = 0; i_p < pPts.size(); ++i_p) {
			for (std::size_t i_T = 0; i_T < TPts.size(); ++i_T) {
				auto pMasses = utils::calculateMasses(m_particleType, m_nf, TPts[i_T]);
				double mu_sq = pMasses.mu*pMasses.mu;
				double mg_sq = pMasses.mg*pMasses.mg;
				double M_sq = pMasses.M*pMasses.M;
				double e = std::sqrt(pPts[i_p]*pPts[i_p] + M_sq);
				double alpha1 = m_alpha_prefactor/std::log(e*TPts[i_T]*utils::constants::INV_LAMBDA_QCD_SQ);

				for (std::size_t i_x = 0; i_x < xPts.size(); ++i_x) {
					m_LdndxTbl[i_tau][i_p][i_T][i_x] = Ldndx(tauPts[i_tau], TPts[i_T], xPts[i_x], mu_sq, mg_sq, M_sq, e, alpha1);
				}
				
				xIntegLimitLow = pMasses.mu*utils::constants::INV_SQRT_2/(pPts[i_p] + e);
				if (m_particleType == utils::ParticleType::Gluon) {
					xIntegLimitHigh = 0.5;
				} else {
					xIntegLimitHigh = 1.0 - pMasses.M/(e + pPts[i_p]);
				}

				m_LNormTbl[i_tau][i_p][i_T] = poly::cubicIntegrate(xPts, m_LdndxTbl[i_tau][i_p][i_T], xIntegLimitLow, xIntegLimitHigh);
			}
		}
	}
}

void LTables::LCollHSeqInit() {
	m_LCollHSeq1.resize(m_LCollMaxPoints);
    m_LCollHSeq2.resize(m_LCollMaxPoints);
    m_LCollHSeq3.resize(m_LCollMaxPoints);
	for (std::size_t i = 0; i < m_LCollMaxPoints; ++i) {
		m_LCollHSeq1[i] = utils::haltonSequence((i+1)*409, 2);
		m_LCollHSeq2[i] = utils::haltonSequence((i+1)*409, 3);
		m_LCollHSeq3[i] = utils::haltonSequence((i+1)*409, 5);
	}
}

std::complex<double> LTables::deltaL2(double q, double w, double mu_sq, double mu_bq) const noexcept{
	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));
	double q_sq = q*q;
	double w_sq = w*w;

	std::complex<double> fn = q_sq + mu_sq*(1.0 - w/2.0/q*log_c);
	fn  = fn*fn;
	fn += (utils::constants::PI_SQ*mu_bq/4.0*w_sq/q_sq);

	return (1.0/fn);
}

std::complex<double> LTables::deltaT2(double q, double w, double mu_sq, double mu_bq) const noexcept {
	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));
	double q_sq = q*q;
	double q_cb = q_sq*q;
	double q_bq = q_cb*q;
	double w_sq = w*w;
	double q_sq_minus_w_sq = q_sq - w_sq;

	std::complex<double> fn = w_sq/q_sq + w*q_sq_minus_w_sq/2.0/q_cb*log_c;
	fn *= (mu_sq/2.0);
	fn += q_sq_minus_w_sq;
	fn = fn*fn;
	fn += (utils::constants::PI_SQ*mu_bq/4.0*w_sq/q_sq*q_sq_minus_w_sq*q_sq_minus_w_sq/4.0/q_bq);

	return (1.0/fn);
}

double LTables::ENumFinite(double p, double T) const noexcept {
	auto pMasses = utils::calculateMasses(m_particleType, m_nf, T);
	double mu_sq = pMasses.mu*pMasses.mu;
	double mu_bq = mu_sq*mu_sq;
	double e = std::sqrt(p*p + pMasses.M * pMasses.M);
	double v = p/e;
	double v_sq = v*v;
	double alpha1 = 4.0*utils::constants::PI/(11.0 - 2.0/3.0*m_nf)/std::log(e*T*utils::constants::INV_LAMBDA_QCD_SQ);
	double alpha2 = 2.0*utils::constants::PI/(11.0 - 2.0/m_nf*3.0)/std::log(pMasses.mu*utils::constants::INV_LAMBDA_QCD);

	//ENumFinite1 integral:
	double ENumFiniteSum1 = 0.0;

	double nfCol;

	double k;
	double kl = 0.0001;
	double kh = m_kmaxColl;
	double kq = kh - kl;
	double expKT;

	double ql = 0.0001;
	double qh, qq, q, q_sq, q_bq, qh1, qh2;
	double qmaxCol = std::sqrt(6.0*e*T);
	
	double wl, wh, wq, w, w_sq;

	#pragma omp parallel for reduction(+:ENumFiniteSum1) private(k,expKT,nfCol, qh,qq,q,q_sq,q_bq, wl,wh,wq,w,w_sq)
	for (std::size_t i=0; i<m_LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + m_LCollHSeq1[i]*kq;
		expKT = std::exp(k/T);
		nfCol = m_Ng/(expKT - 1.0) + m_nf/(expKT + 1.0);

		qh = ((qmaxCol < k) ? qmaxCol : k);
		qq = qh - ql;
		q = ql + m_LCollHSeq2[i]*qq;
		q_sq = q*q;
		q_bq = q_sq*q_sq;

		wl = -q;
		wh = q;
		wq = wh - wl;
		w = wl + m_LCollHSeq3[i]*wq;
		w_sq = w*w;

		fn_comp  = 2.0*utils::constants::INV_HBARC_GEVFM*m_CR*alpha1*alpha2*utils::constants::INV_PI/v_sq*nfCol*w*utils::unitStep(v_sq*q_sq - w_sq);
		fn_comp *= (deltaL2(q, w, mu_sq, mu_bq)*((2.0*k + w)*(2.0*k + w) - q_sq)/2.0 + deltaT2(q, w, mu_sq, mu_bq)*(q_sq - w_sq)/4.0/q_bq*((2.0*k + w)*(2.0*k + w) + q_sq)*(v_sq*q_sq - w_sq));

		ENumFiniteSum1 += fn_comp.real()*qq*wq;
	}

	ENumFiniteSum1 = ENumFiniteSum1*kq/static_cast<double>(m_LCollMaxPoints);
	
	//ENumFinite2 integral:
	double ENumFiniteSum2 = 0.0;

	#pragma omp parallel for reduction(+:ENumFiniteSum2) private(k,expKT,nfCol, ql,qh,qh1,qh2,qq,q,q_sq,q_bq, wl,wh,wq,w,w_sq)
	for (std::size_t i=0; i<m_LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + m_LCollHSeq1[i]*kq;
		expKT = std::exp(k/T);
		nfCol = m_Ng/(expKT - 1.0) + m_nf/(expKT + 1.0);

		ql = ((qmaxCol < k) ? qmaxCol : k);
		qh1 = 2.0*k*(1.0 - k/e)/(1.0 - v + 2.0*k/e);
		qh2 = ((k > qh1) ? k : qh1);
		qh = ((qmaxCol < qh2) ? qmaxCol : qh2);
		qq = qh - ql;
		q = ql + m_LCollHSeq2[i]*qq;
		q_sq = q*q;
		q_bq = q_sq*q_sq;

		wl = q - 2.0*k;
		wh = q;
		wq = wh - wl;
		w = wl + m_LCollHSeq3[i] * wq;
		w_sq = w*w;

		fn_comp  = 2.0*utils::constants::INV_HBARC_GEVFM*m_CR*alpha1*alpha2*utils::constants::INV_PI/v_sq*nfCol*w*utils::unitStep(v_sq*q_sq - w_sq);
		fn_comp *= (deltaL2(q, w, mu_sq, mu_bq)*((2.0*k + w)*(2.0*k + w) - q_sq)/2.0 + deltaT2(q, w, mu_sq, mu_bq)*(q_sq - w_sq)/4.0/q_bq*((2.0*k + w)*(2.0*k + w) + q_sq)*(v_sq*q_sq - w_sq));

		ENumFiniteSum2 += fn_comp.real()*qq*wq;
	}
	
	ENumFiniteSum2 = ENumFiniteSum2*kq/static_cast<double>(m_LCollMaxPoints);

	return (ENumFiniteSum1 + ENumFiniteSum2);
}

void LTables::CollLTables()
{
	LCollHSeqInit();

	const std::vector<double>& pCollPts = m_Grids.pCollPts();
	const std::vector<double>& TCollPts = m_Grids.TCollPts();

	m_LCollTbl.resize(
		pCollPts.size(), std::vector<double>(
			TCollPts.size(), 0.0
		)
	);

	for (std::size_t i_p = 0; i_p < pCollPts.size(); ++i_p) {
		for (std::size_t i_T = 0; i_T < TCollPts.size(); ++i_T) {
			m_LCollTbl[i_p][i_T] = ENumFinite(pCollPts[i_p], TCollPts[i_T]);
		}
	}
}

int LTables::exportLTables() const
{
	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	{//exporting Ldndx table
		const std::string path_out = "./ltables/ldndx_nf=" + nfss.str() + "_" + m_pName + "_xB=" + xBss.str() + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open Ldndx table export file." << std::endl;
			return -1;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(12) <<   "tau" << " ";
		file_out << std::fixed << std::setw(14) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(12) <<     "x" << " ";
		file_out << std::fixed << std::setw(17) << "Ldndx" << "\n";

		for (std::size_t itau=0; itau<m_Grids.tauPtsLength(); itau++) {
			for (std::size_t ip=0; ip<m_Grids.pPtsLength(); ip++) {
				for (std::size_t iT=0; iT<m_Grids.TPtsLength(); iT++) {
					for (std::size_t ix=0; ix<m_Grids.xPtsLength(); ix++) {
						file_out << std::fixed 		<< std::setw(13) << std::setprecision(10) << m_Grids.tauPts(itau) << " ";
						file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pPts(ip) << " ";
						file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TPts(iT) << " ";
						file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.xPts(ix) << " ";
						file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LdndxTbl[itau][ip][iT][ix] << "\n";
					}
				}
			}
		}

		file_out.close();
	}

	{//exporting LNorm table
		const std::string path_out = "./ltables/lnorm_nf=" + nfss.str() + "_" + m_pName + "_xB=" + xBss.str() + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open LNorm table export file." << std::endl;
			return -2;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(12) <<   "tau" << " ";
		file_out << std::fixed << std::setw(14) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(17) << "LNorm" << "\n";

		for (std::size_t itau=0; itau<m_Grids.tauPtsLength(); itau++) {
			for (std::size_t ip=0; ip<m_Grids.pPtsLength(); ip++) {
				for (std::size_t iT=0; iT<m_Grids.TPtsLength(); iT++) {
					file_out << std::fixed 		<< std::setw(13) << std::setprecision(10) << m_Grids.tauPts(itau) << " ";
					file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pPts(ip) << " ";
					file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TPts(iT) << " ";
					file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LNormTbl[itau][ip][iT] << "\n";
				}
			}
		}

		file_out.close();
	}

	{//exporting LColl table
		std::string path_out = "./ltables/lcoll_nf=" + nfss.str() + "_" + m_pName + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open LColl table export file." << std::endl;
			return -3;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(13) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(17) << "LColl" << "\n";

		for (std::size_t ip=0; ip<m_Grids.pCollPtsLength(); ip++) {
			for (std::size_t iT=0; iT<m_Grids.TCollPtsLength(); iT++) {
				file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pCollPts(ip) << " ";
				file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TCollPts(iT) << " ";
				file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LCollTbl[ip][iT] << "\n";
			}
		}

		file_out.close();
	}

	return 1;
}