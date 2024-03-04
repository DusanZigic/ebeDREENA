#include "mainheader.hpp"
#include "ltables.hpp"
#include "grids.hpp"
#include "polyintegration.hpp"
#include "arsenal.hpp"
#include "importexport.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

std::string LT_sNN;  	//collision energy
double LT_nf;		    //effective number of flavours
std::string LT_pName;	//particle name
double LT_xB;		    //xB value
gridPoints LT_Grids;    //grids
size_t LdndxMaxPoints;  //maximal number of points for Ldndx integration
size_t LCollMaxPoints;  //maximal number of points for collisional integration
double LT_TCRIT;		//critical temperature

static const double Ng = 3.0;		//effective number of gluons
static const double lambda = 0.2;   //QCD scale
static const double kmaxColl = 5.0; //kMaxColl value
static       double LT_CR;		    //Casimir (3 for gluons, 4/3 for quakrs)


static std::vector<double> LdndxHSeq1, LdndxHSeq2, LdndxHSeq3;
static void LdndxHSeqInit()
{
	for (size_t i=0; i<LdndxMaxPoints; i++) {
		LdndxHSeq1.push_back(haltonSequence((i+1)*409, 2));
		LdndxHSeq2.push_back(haltonSequence((i+1)*409, 3));
		LdndxHSeq3.push_back(haltonSequence((i+1)*409, 5));
	}
}

static double dElossDYN(double tau, double p, double x, double k, double q, double varphi, double T)
{
	double mu = 0.197*std::sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0 + LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));
	double mg = mu / std::sqrt(2.0);
	double M = 0.0;
	if (LT_pName == "Bottom") M = 4.75;
	else if (LT_pName == "Charm") M = 1.2;
	else if (LT_pName == "Gluon") M = mu/std::sqrt(2.0);
	else M = mu/std::sqrt(6.0);

	double b = std::sqrt(mg*mg + M * M*x*x);
	double e = std::sqrt(p*p + M * M);
	double alpha  = 4.0*M_PI/(11.0 - 2.0*LT_nf/3.0)/std::log((k*k + mg*mg + M*M*x*x)/x/lambda/lambda);
	double alpha1 = 4.0*M_PI/(11.0 - 2.0*LT_nf/3.0)/std::log(e*T/0.2/0.2);

	double fn = 1.0;
	fn *= 1.0 / 0.197*LT_CR*alpha/M_PI*3.0*alpha1*T*2.0*k*q/M_PI;
	fn *= (mu*mu - mu*mu*LT_xB*LT_xB)/(q*q + mu*mu*LT_xB*LT_xB)/(q*q + mu*mu);

	double psi = (k*k + q*q + 2.0*k*q*std::cos(varphi) + b*b)/2.0/x/e*tau/0.197;

	fn *= (1 - std::cos(psi));
	fn *= 2.0/(k*k + b*b)/(k*k + q*q + 2.0*k*q*cos(varphi) + b*b)/(k*k + q*q + 2.0*k*q*std::cos(varphi) + b*b);
	fn *= (-1.0*k*q*std::cos(varphi)*(k*k + q*q + 2.0*k*q*std::cos(varphi)) + b*b*(k*q*std::cos(varphi) + q*q));

	return fn;
}

static double Ldndx(double tau, double p, double T, double x)
{
	double mu = 0.197*std::sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));
	double mg = mu / std::sqrt(2.0);
	double M = 0.0;
	if (LT_pName == "Bottom") M = 4.75;
	else if (LT_pName == "Charm") M = 1.2;
	else if (LT_pName == "Gluon") M = mg;
	else M = mu/sqrt(6.0);
	double e = sqrt(p*p + M * M);

	double kl = 0.00000001; 
	double kh = 2.0*x*(1 - x)*e;
	double kq = (kh - kl);
	double ql = 0.000001;
	double qh = sqrt(4.0*e*T);
	double qq = qh - ql;
	double phil = 0.0;
	double phih = M_PI;
	double phiq = (phih - phil);
	double sum = 0.0; //integration sum
	double k, q, phi; //integration variables

	#pragma omp parallel for reduction(+:sum) private(k,q,phi)
	for (size_t i = 0; i<LdndxMaxPoints; i++) {
		  k  =   kl + LdndxHSeq1[i]*kq;
		  q  =   ql + LdndxHSeq2[i]*qq;
		phi  = phil + LdndxHSeq3[i]*phiq;
		sum += 2*dElossDYN(tau, p, x, k, q, phi, T)/x;
	}

	return (sum*kq*qq*phiq/LdndxMaxPoints);
}

static void RadLTables(std::vector<std::vector<std::vector<std::vector<double>>>> &ldndxtbl, std::vector<std::vector<std::vector<double>>> &lnormtbl)
{
	LdndxHSeqInit();

	double tau, p, T, x, mu, M, xIntegLimitLow, xIntegLimitHigh;

	for (size_t itau=0; itau<LT_Grids.tauPtsLength(); itau++) {
		tau = LT_Grids.tauPts(itau);

		for (size_t ip=0; ip<LT_Grids.pPtsLength(); ip++) {
			p = LT_Grids.pPts(ip);

			for (size_t iT=0; iT<LT_Grids.TPtsLength(); iT++) {
				T = LT_Grids.TPts(iT);
				
				for (size_t ix=0; ix<LT_Grids.xPtsLength(); ix++) {
					x = LT_Grids.xPts(ix);
					ldndxtbl[itau][ip][iT][ix] = Ldndx(tau, p, T, x);
				}
				
				mu = 0.197*std::sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));
				if (LT_pName == "Bottom") M = 4.75;
				else if (LT_pName == "Charm") M = 1.2;
				else if (LT_pName == "Gluon") M = mu/sqrt(2.0);
				else M = mu/sqrt(6.0);
				
				xIntegLimitLow = mu/std::sqrt(2.0)/(p + std::sqrt(p*p + M*M));
				if (LT_pName == "Gluon") xIntegLimitHigh = 0.5;
				else xIntegLimitHigh = 1.0 - M/(std::sqrt(p*p + M*M) + p);

				lnormtbl[itau][ip][iT] = cubicIntegrate(LT_Grids.xPts(), ldndxtbl[itau][ip][iT], xIntegLimitLow, xIntegLimitHigh);
			}
		}
	}
}

static std::vector<double> LCollHSeq1, LCollHSeq2, LCollHSeq3;
static void LCollHSeqInit()
{
	for (size_t i=0; i<LCollMaxPoints; i++) {
		LCollHSeq1.push_back(haltonSequence((i+1)*409, 2));
		LCollHSeq2.push_back(haltonSequence((i+1)*409, 3));
		LCollHSeq3.push_back(haltonSequence((i+1)*409, 5));
	}
}

static std::complex<double> deltaL2(double q, double w, double T)
{
	double mu = 0.197*sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));

	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));

	std::complex<double> fn = q*q + mu*mu*(1.0 - w/2.0/q*log_c);
	fn  = fn*fn;
	fn += (M_PI*M_PI*mu*mu*mu*mu/4.0*w*w/q/q);

	return (1.0/fn);
}

static std::complex<double> deltaT2(double q, double w, double T)
{
	double mu = 0.197*sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));

	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));

	std::complex<double> fn = w*w/q/q + w*(q*q - w*w)/2.0/q/q/q*log_c;
	fn *= (mu*mu/2.0);
	fn += (q*q - w*w);
	fn = fn*fn;
	fn += (M_PI*M_PI*mu*mu*mu*mu/4.0*w*w/q/q*(q*q - w*w)*(q*q - w*w)/4.0/q/q/q/q);

	return (1.0/fn);
}

static double ENumFinite(double p, double T)
{
	double mu = 0.197*sqrt((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda/productLog((-8.0*(6.0+LT_nf)*M_PI*M_PI*T*T)/(2.0*LT_nf-33.0)/lambda/lambda));
	double M = 1.0;
	if (LT_pName == "Bottom") M = 4.75;
	else if (LT_pName == "Charm") M = 1.2;
	else if (LT_pName == "Gluon") M = mu/std::sqrt(2.0);
	else M = mu/std::sqrt(6.0);
	double e = std::sqrt(p*p + M * M);
	double v = p/e;
	double alpha1 = 4.0*M_PI/(11.0 - 2.0/3.0*LT_nf)/std::log(e*T/0.2/0.2);
	double alpha2 = 2.0*M_PI/(11.0 - 2.0/LT_nf*3.0)/std::log(mu/0.2);

	//ENumFinite1 integral:
	double ENumFiniteSum1 = 0.0;

	double nfCol;

	double k;
	double kl = 0.0001;
	double kh = kmaxColl;
	double kq = kh - kl;

	double ql = 0.0001;
	double qh, qq, q, qmaxCol, qh1, qh2;

	double wl, wh, wq, w;

	#pragma omp parallel for reduction(+:ENumFiniteSum1) private(k,nfCol, qh,qq,q,qmaxCol, wl,wh,wq,w)
	for (size_t i=0; i<LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + LCollHSeq1[i]*kq;
		nfCol = Ng/(std::exp(k/T) - 1.0) + LT_nf/(std::exp(k/T) + 1.0);

		qmaxCol = std::sqrt(6.0*e*T);
		qh = ((qmaxCol < k) ? qmaxCol : k);
		qq = qh - ql;
		q = ql + LCollHSeq2[i]*qq;

		wl = -q;
		wh = q;
		wq = wh - wl;
		w = wl + LCollHSeq3[i]*wq;

		fn_comp  = 2.0/0.197*LT_CR*alpha1*alpha2/M_PI/v/v*nfCol*w*unitStep(v*v*q*q - w*w);
		fn_comp *= (deltaL2(q, w, T)*((2.0*k + w)*(2.0*k + w) - q*q)/2.0 + deltaT2(q, w, T)*(q*q - w*w)/4.0/q/q/q/q*((2.0*k + w)*(2.0*k + w) + q*q)*(v*v*q*q - w*w));

		ENumFiniteSum1 += fn_comp.real()*qq*wq;
	}

	ENumFiniteSum1 = ENumFiniteSum1*kq/LCollMaxPoints;

	
	//ENumFinite2 integral:
	double ENumFiniteSum2 = 0.0;

	#pragma omp parallel for reduction(+:ENumFiniteSum2) private(k,nfCol, ql,qh,qh1,qh2,qq,q,qmaxCol, wl,wh,wq,w)
	for (size_t i=0; i<LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + LCollHSeq1[i]*kq;
		nfCol = Ng/(std::exp(k/T) - 1.0) + LT_nf/(std::exp(k/T) + 1.0);

		qmaxCol = std::sqrt(6.0*e*T);
		ql = ((qmaxCol < k) ? qmaxCol : k);
		qh1 = 2.0*k*(1.0 - k/e)/(1.0 - v + 2.0*k/e);
		qh2 = ((k > qh1) ? k : qh1);
		qh = ((qmaxCol < qh2) ? qmaxCol : qh2);
		qq = qh - ql;
		q = ql + LCollHSeq2[i]*qq;

		wl = q - 2.0*k;
		wh = q;
		wq = wh - wl;
		w = wl + LCollHSeq3[i] * wq;

		fn_comp  = 2.0/0.197*LT_CR*alpha1*alpha2/M_PI/v/v*nfCol*w*unitStep(v*v*q*q - w*w);
		fn_comp *= (deltaL2(q, w, T)*((2.0*k + w)*(2.0*k + w) - q*q)/2.0 + deltaT2(q, w, T)*(q*q - w*w)/4.0/q/q/q/q*((2.0*k + w)*(2.0*k + w) + q*q)*(v*v*q*q - w*w));

		ENumFiniteSum2 += fn_comp.real()*qq*wq;
	}
	
	ENumFiniteSum2 = ENumFiniteSum2*kq/LCollMaxPoints;


	return (ENumFiniteSum1 + ENumFiniteSum2);
}

static void CollLTables(std::vector<std::vector<double>> &lcolltbl)
{
	LCollHSeqInit();

	for (size_t ip=0; ip<LT_Grids.pCollPtsLength(); ip++) {
		for (size_t iT=0; iT<LT_Grids.TCollPtsLength(); iT++) {
			lcolltbl[ip][iT] = ENumFinite(LT_Grids.pCollPts(ip), LT_Grids.TCollPts(iT));
		}
	}
}

void GenerateLTables()
{
	LT_nf = 3.0; if (LT_sNN == "200GeV") LT_nf = 2.5;

	LT_Grids.setGridPoints(LT_sNN, LT_pName, LT_TCRIT);

    LT_CR = LT_pName == "Gluon" ? 3.0 : 4.0/3.0;

	std::vector<std::vector<std::vector<std::vector<double>>>> LdndxTbl(LT_Grids.tauPtsLength(),
																		std::vector<std::vector<std::vector<double>>>(LT_Grids.pPtsLength(),
																		std::vector<std::vector<double>>(LT_Grids.TPtsLength(),
																		std::vector<double>(LT_Grids.xPtsLength(), 0.0))));;
	std::vector<std::vector<std::vector<double>>> LNormTbl(LT_Grids.tauPtsLength(),
														   std::vector<std::vector<double>>(LT_Grids.pPtsLength(),
														   std::vector<double>(LT_Grids.TPtsLength(), 0.0)));
    RadLTables(LdndxTbl, LNormTbl);

	std::vector<std::vector<double>> LCollTbl(LT_Grids.pCollPtsLength(), std::vector<double>(LT_Grids.TCollPtsLength(), 0.0));
    CollLTables(LCollTbl);

	if (exportLTables(LdndxTbl, LNormTbl, LCollTbl) != 1) return;
}