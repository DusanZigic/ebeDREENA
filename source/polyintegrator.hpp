#ifndef POLY_INTINTEGRATOR_HPP
#define POLY_INTINTEGRATOR_HPP

#include <vector>
#include <cmath>

/**
 * SANITY CHECK STRATEGY:
 * By default, this library is optimized for "hot-path" execution (no checks, maximum speed).
 * This is the high-performance path, allowing the compiler to use branchless 
 * instructions and vectorization.
 * 
 * If -DFAST_SIM_STRICT is defined, the following safety gates are enabled:
 * 1. VALIDATE_DATA_SIZE: ensures the number of data points is larger or equal to 2 for linear
 *    and 4 for cubic integration.
 * 2. VALIDATE_POLY_RANGE: ensures low and high integration parameters are within x data range.
 */
#ifdef FAST_SIM_STRICT
    #include <iostream>
    #include <stdexcept>
    #include <string>

    #define VALIDATE_DATA_SIZE(objName, size, required)                                                      \
        if ((size) < (required)) {                                                                           \
            std::string msg = "\n[DATA SIZE ERROR] " + std::string(objName) +                                \
                            " | Provided: " + std::to_string((size)) +                                       \
                            " | Required: " + std::to_string((required));                                    \
            std::cerr << msg << std::endl;                                                                   \
            throw std::invalid_argument(msg);                                                                \
        }

    #define VALIDATE_POLY_RANGE(objName, val, minV, maxV)                                                    \
        if ((val) < (minV) || (val) > (maxV)) {                                                              \
            std::string msg = "\n[RANGE ERROR] " + std::string((objName)) +                                  \
                              " | Value: " + std::to_string((val)) +                                         \
                              " | Limits: [" + std::to_string((minV)) + ", " + std::to_string((maxV)) + "]"; \
            std::cerr << msg << std::endl;                                                                   \
            throw std::out_of_range(msg);                                                                    \
        }
#else
    #define VALIDATE_POLY_RANGE(objName, val, minV, maxV)
    #define VALIDATE_DATA_SIZE(objName, size, required)
#endif

namespace poly {

namespace detail {
    template <typename T>
    inline std::size_t locatePoint(const std::vector<T>& data, T x, int order) {
        int n = static_cast<int>(data.size());
        
        int jl = 0;
        int ju = n - 1;
        bool ascnd = (data.back() >= data.front());

        while (ju - jl > 1) {
            int jm = (ju + jl) >> 1;
            if ((x >= data[jm]) == ascnd) jl = jm;
            else ju = jm;
        }
        
        int mm = order + 1;
        int pointLocation = std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
        return static_cast<std::size_t>(pointLocation);
    }

    template <typename T>
    inline void getCubicCoeffs(const T* x, const T* f, double* coeff) {
        T s[4] = {static_cast<T>(0)};;
        T phi, ff, b;
        s[3] = -x[0];

        for (int i = 1; i < 4; i++) {
            for (int j = 4 - 1 - i; j < 4 - 1; j++)
                s[j] -= x[i] * s[j+1];
            s[3] -= x[i];
        }

        for (int j = 0; j < 4; j++) {
            phi = 4.0;
            for (int k = 3; k > 0; k--)
                phi = k * s[k] + x[j] * phi;
            
            ff = f[j] / phi;        
            b = 1.0;
            for (int k = 3; k >= 0; k--) {
                coeff[k] += b * ff;
                b = s[k] + x[j] * b;
            }
        }
    }

} // namespace detail

template <typename T>
inline T linearIntegrate(const std::vector<T> &x, const std::vector<T> &f) {
    VALIDATE_DATA_SIZE("LinearIntegrate", x.size(), 2);
    
    T total = static_cast<T>(0);
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        T xL = x[i];
        T xR = x[i + 1];
        T dx = xR - xL;
        T k = (f[i + 1] - f[i]) / dx;
        T c = f[i] - k * xL;

        total += static_cast<T>(0.5) * k * (xR * xR - xL * xL) + c * (xR - xL);
    }
    return total;
}

template <typename T>
inline T cubicIntegrate(const std::vector<T> &x, const std::vector<T> &f) {
    VALIDATE_DATA_SIZE("CubicIntegrate", x.size(), 4);
    
    T total = static_cast<T>(0);
    for (std::size_t i = 0; i < x.size() - 1; ++i) {
        std::size_t winIdx = detail::locatePoint(x, x[i], 3);
        T coeffs[4] = {static_cast<T>(0)};
        detail::getCubicCoeffs(&x[winIdx], &f[winIdx], coeffs);

        T xL = x[i];
        T xR = x[i + 1];
        T pL = xL; 
        T pR = xR;
        
        for (int j = 0; j < 4; ++j) {
            total += (coeffs[j] / static_cast<T>(j + 1)) * (pR - pL);
            pL *= xL;
            pR *= xR;
        }
    }
    return total;
}

template <typename T>
inline T linearIntegrate(const std::vector<T> &x, const std::vector<T> &f, T low, T high) {
    VALIDATE_DATA_SIZE("LinearIntegrate", xdata.size(), 2);
    VALIDATE_POLY_RANGE("LinearIntegrate", low, x.front(), x.back());
    VALIDATE_POLY_RANGE("LinearIntegrate", high, x.front(), x.back());

    std::size_t iStart = detail::locatePoint(x, low, 1);
    std::size_t iEnd = detail::locatePoint(x, high, 1);
    T total = static_cast<T>(0);;

    for (std::size_t i = iStart; i <= iEnd; ++i) {
        T xL = (i == iStart) ? low : x[i];
        T xR = (i == iEnd) ? high : x[i+1];
        
        T dx = x[i+1] - x[i];
        T k = f[i+1] - f[i] / dx;
        T c = f[i] - k * x[i];

        total += static_cast<T>(0.5) * k * (xR * xR - xL * xL) + c * (xR - xL);
    }
    return total;
}

template <typename T>
inline T cubicIntegrate(const std::vector<T> &x, const std::vector<T> &f, T low, T high) {
    VALIDATE_DATA_SIZE("CubicIntegrate", xdata.size(), 4);
    VALIDATE_POLY_RANGE("CubicIntegrate", low, x.front(), x.back());
    VALIDATE_POLY_RANGE("CubicIntegrate", high, x.front(), x.back());

    std::size_t iStart = detail::locatePoint(x, low, 1);
    std::size_t iEnd = detail::locatePoint(x, high, 1);
    T total = static_cast<T>(0);

    for (std::size_t i = iStart; i <= iEnd; ++i) {
        T xL = (i == iStart) ? low : x[i];
        T xR = (i == iEnd) ? high : x[i+1];

        std::size_t winIdx = detail::locatePoint(x, x[i], 3);
        T coeffs[4] = {static_cast<T>(0)};
        detail::getCubicCoeffs(&x[winIdx], &f[winIdx], coeffs);

        T pL = xL; 
        T pR = xR;
        for (int j = 0; j < 4; ++j) {
            total += (coeffs[j] / static_cast<T>(j + 1)) * (pR - pL);
            pL *= xL;
            pR *= xR;
        }
    }
    return total;
}

} // namespace poly

#endif //POLY_INTINTEGRATOR_HPP