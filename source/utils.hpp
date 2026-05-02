#include <cmath>

namespace {
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
} // namespace

namespace utils {
    const double HBARC_GEVFM = 0.197;

    inline double debyeMass(double nf, double lambda, double T) {
        double numerator = -8.0*(6.0 + nf)*M_PI*M_PI*T*T;
        double denominator = (2.0*nf - 33.0)*lambda*lambda;
        double x = numerator / denominator;
        double mu_squared  = x / productLog(x);
        return HBARC_GEVFM * std::sqrt(mu_squared);
    }

    inline double unitStep(double x) {
        return (x < 0.0) ? 0.0 : 1.0;
    }

    inline long double unitStep(long double x) {
        return (x < 0.0L) ? 0.0L : 1.0L;
    }

} // namspace utils