#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <cmath>

namespace utils {

    namespace constants {
        constexpr double HBARC_GEVFM     = 0.197;             // reduced Planck's constant (hbar) times speed of light (c) in units of GeV·fm
        constexpr double INV_HBARC_GEVFM = 1.0 / HBARC_GEVFM; // precalculated inverse HBARC_GEVFM

        constexpr double LAMBDA_QCD        = 0.2;                     // QCD scale (lambda)
        constexpr double LAMBDA_QCD_SQ     = LAMBDA_QCD * LAMBDA_QCD; // precalculated squared lambda
        constexpr double INV_LAMBDA_QCD    = 1.0 / LAMBDA_QCD;        // precalculated inverse QCD scale (lambda)
        constexpr double INV_LAMBDA_QCD_SQ = 1.0 / LAMBDA_QCD_SQ;     // precalculated inverse squared lambda

        constexpr double M_BOTTOM = 4.75; // bottom quark mass
        constexpr double M_CHARM  = 1.20; // charm quark mass

        constexpr double PI     = 3.14159265358979323846; // pi
        constexpr double PI_SQ  = PI * PI;                // precalculated pi squared
        constexpr double INV_PI = 1.0 / PI;               // precalculated inverse pi

        constexpr double INV_SQRT_2 = 0.70710678118654746172; // precalculated 1 / sqrt(2)
        constexpr double INV_SQRT_6 = 0.40824829046386307274; // precalculated 1 / sqrt(6)

    } //namespace constants

    enum class ParticleType: std::uint8_t { Bottom, Charm, Gluon, LQuarks }; // particle type enum class
    struct ParticleMasses { double mu; double mg; double M; };               // masses structure containing Debye, gluon and jet mass

    namespace detail {
        inline double productLog(double x) {
            if (x == 0.0) {
    		    return 0.0;
            }

            double w0, w1;
            if (x > 0.0) {
                w0 = std::log(1.2 * x / std::log(2.4 * x / std::log1p(2.4 * x)));
            }
            else {
                double v = 1.4142135623730950488 * std::sqrt(1.0 + 2.7182818284590452354 * x);
                double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
                double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
                w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
            }

            while (true) {
                double e = std::exp(w0);
                double f = w0 * e - x;
                w1 = w0 - f / ((e * (w0 + 1.0) - (w0 + 2.0) * f / (w0 + w0 + 2.0)));
                if (std::abs(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
                    break;
                }
                w0 = w1;
            }
            return w1;
        }
    } //namespace detail

    inline double debyeMass(double nf, double T) {
        double numerator = -8.0*(6.0 + nf)*constants::PI_SQ*T*T;
        double denominator = (2.0*nf - 33.0)*constants::LAMBDA_QCD_SQ;
        double x = numerator / denominator;
        double mu_squared  = x / detail::productLog(x);
        return constants::HBARC_GEVFM * std::sqrt(mu_squared);
    }

    inline ParticleMasses calculateMasses(ParticleType particleType, double nf, double T) {
        double mu = utils::debyeMass(nf, T);
        double mg = mu * constants::INV_SQRT_2;
        double M = 0.0;

        switch (particleType) {
            case ParticleType::Bottom:  M = constants::M_BOTTOM;      break;
            case ParticleType::Charm:   M = constants::M_CHARM;       break;
            case ParticleType::Gluon:   M = mg;                       break;
            case ParticleType::LQuarks: M = mu*constants::INV_SQRT_6; break;
        }

        return {mu, mg, M};
    }

    inline double unitStep(double x) {
        return (x < 0.0) ? 0.0 : 1.0;
    }

    inline long double unitStep(long double x) {
        return (x < 0.0L) ? 0.0L : 1.0L;
    }

    inline double haltonSequence(int index, int base) {
        double f = 1.0, res = 0.0;
        while (index > 0) {
            f = f / static_cast<double>(base);
            res += f * static_cast<double>(index % base);
            index = index / base; // integer division
        }
        return res;
    }

} // namspace utils

#endif // UTILS_HPP