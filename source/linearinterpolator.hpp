#ifndef LINEAR_INTERPOLATOR_HPP
#define LINEAR_INTERPOLATOR_HPP

#include <vector>
#include <array>
#include <cstddef>

/**
 * SANITY CHECK STRATEGY:
 * By default, this library is optimized for "hot-path" execution (no checks, maximum speed).
 * This is the high-performance path, allowing the compiler to use branchless 
 * instructions and vectorization.
 * 
 * If -DFAST_SIM_STRICT is defined, the following safety gates are enabled:
 * 1. VALIDATE_DIMENSIONS: ensures the number of arguments in interpolate() 
 *    matches the table dimensions (prevents logic errors).
 * 2. VALIDATE_BOUNDS: throws an exception on extrapolation
 */
#ifdef FAST_SIM_STRICT
    #include <iostream>
    #include <stdexcept>
    #include <string>

    #define VALIDATE_DIMENSIONS(expected, actual)                                                        \
        if ((expected) != (actual)) {                                                                    \
            std::string msg = "\n[DIMENSION ERROR]"                                                      \
                              " | Object: " + std::string(m_name) +                                      \
                              " | Expected: " + std::to_string(expected) + "D" +                         \
                              " | Received: " + std::to_string(actual) + "D";                            \
            std::cerr << msg << std::endl;                                                               \
            throw std::invalid_argument(msg);                                                            \
        }

    #define VALIDATE_BOUNDS(d, val, minV, maxV)                                                          \
        if ((val) < (minV) || (val) > (maxV)) {                                                          \
            std::string msg = "\n[EXTRAPOLATION ERROR]"                                                  \
                              " | Object: " + std::string(m_name) +                                      \
                              " | Dimension: " + std::to_string(d) +                                     \
                              " | Value: " + std::to_string(val) +                                       \
                              " | Domain: [" + std::to_string(minV) + ", " + std::to_string(maxV) + "]"; \
            std::cerr << msg << std::endl;                                                               \
			throw std::out_of_range(msg);                                                                \
        }
#else
    #define VALIDATE_DIMENSIONS(expected, actual)
    #define VALIDATE_BOUNDS(d, val, minV, maxV) 
#endif

template<typename T>
class LinearInterpolator {
public:
	LinearInterpolator(const char* name = "unnamed") : m_name(name), m_dim(0) {}

	LinearInterpolator(const std::vector<T> &x, const std::vector<T> &f, const char* name = "unnamed") : m_name(name) {
		initializeInternal({x}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f, const char* name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<std::vector<T>> &f_grid, const char* name = "unnamed") : m_name(name) {
		std::vector<T> f_flat;
        f_flat.reserve(x1.size() * x2.size());
        for (const auto& row : f_grid) {
            f_flat.insert(f_flat.end(), row.begin(), row.end());
        }
        initializeInternal({x1, x2}, f_flat);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f, const char* name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2, x3}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f, const char* name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2, x3, x4}, f);
	}

	void setData(const std::vector<T> &x, const std::vector<T> &f, const char* name = "unnamed") {
		m_name = name;
		initializeInternal({x}, f);
	}

	void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f, const char* name = "unnamed") {
		m_name = name;
		initializeInternal({x1, x2}, f);
	}

    void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f, const char* name = "unnamed") {
		m_name = name;
		initializeInternal({x1, x2, x3}, f);
	}

    void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f, const char* name = "unnamed") {
		m_name = name;
		initializeInternal({x1, x2, x3, x4}, f);
	}

	inline T interpolate(T x) const {
        VALIDATE_DIMENSIONS(m_dim, 1);

        std::size_t i = locateGridIndex(0, x);

        std::size_t ss = m_searchStrides[0];
        
        T x0 = m_axes[0][i * ss], x1 = m_axes[0][(i + 1) * ss];
        T frac = (x - x0) / (x1 - x0);
        return m_values[i] * (1 - frac) + m_values[i + 1] * frac;
    }

    inline T interpolate(T x1, T x2) const {
        VALIDATE_DIMENSIONS(m_dim, 2);

        std::size_t i1 = locateGridIndex(0, x1), i2 = locateGridIndex(1, x2);

		std::size_t ss1 = m_searchStrides[0], ss2 = m_searchStrides[1];
        std::size_t s0 = m_memStrides[0];

        T w1 = (x1 - m_axes[0][i1 * ss1]) / (m_axes[0][(i1 + 1) * ss1] - m_axes[0][i1 * ss1]);
        T w2 = (x2 - m_axes[1][i2 * ss2]) / (m_axes[1][(i2 + 1) * ss2] - m_axes[1][i2 * ss2]);

        std::size_t idx = i1 * s0 + i2;
        T v00 = m_values[idx];
        T v01 = m_values[idx + 1];
        T v10 = m_values[idx + s0];
        T v11 = m_values[idx + s0 + 1];

        return (1-w1)*(1-w2)*v00 + (1-w1)*w2*v01 + w1*(1-w2)*v10 + w1*w2*v11;
    }

    inline T interpolate(T x1, T x2, T x3) const {
        VALIDATE_DIMENSIONS(m_dim, 3);
        
        std::size_t i1 = locateGridIndex(0, x1), i2 = locateGridIndex(1, x2), i3 = locateGridIndex(2, x3);

        std::size_t ss1 = m_searchStrides[0], ss2 = m_searchStrides[1], ss3 = m_searchStrides[2];
        std::size_t s0 = m_memStrides[0], s1 = m_memStrides[1];

        T w1 = (x1 - m_axes[0][i1*ss1]) / (m_axes[0][(i1+1)*ss1] - m_axes[0][i1*ss1]);
        T w2 = (x2 - m_axes[1][i2*ss2]) / (m_axes[1][(i2+1)*ss2] - m_axes[1][i2*ss2]);
        T w3 = (x3 - m_axes[2][i3*ss3]) / (m_axes[2][(i3+1)*ss3] - m_axes[2][i3*ss3]);

        auto g = [&](std::size_t o1, std::size_t o2, std::size_t o3) { 
            return m_values[(i1+o1)*s0 + (i2+o2)*s1 + (i3+o3)]; 
        };

        T c00 = g(0,0,0)*(1-w1) + g(1,0,0)*w1;
        T c01 = g(0,0,1)*(1-w1) + g(1,0,1)*w1;
        T c10 = g(0,1,0)*(1-w1) + g(1,1,0)*w1;
        T c11 = g(0,1,1)*(1-w1) + g(1,1,1)*w1;

        T c0 = c00*(1-w2) + c10*w2;
        T c1 = c01*(1-w2) + c11*w2;

        return c0*(1-w3) + c1*w3;
    }

    inline T interpolate(T x1, T x2, T x3, T x4) const {
        VALIDATE_DIMENSIONS(m_dim, 4);
        
        std::size_t i[4] = {locateGridIndex(0, x1), locateGridIndex(1, x2), locateGridIndex(2, x3), locateGridIndex(3, x4)};

        T w[4];
        for(int d=0; d<4; ++d) {
            std::size_t base = i[d] * m_searchStrides[d];
            std::size_t next = (i[d] + 1) * m_searchStrides[d];
            w[d] = ( (d==0?x1:(d==1?x2:(d==2?x3:x4))) - m_axes[d][base] ) / (m_axes[d][next] - m_axes[d][base]);
        }

        T res = 0;
        for (int b = 0; b < 16; ++b) {
            T wt = 1.0;
            std::size_t idx = 0;
            for (int d = 0; d < 4; ++d) {
                bool bit = (b >> (3 - d)) & 1;
                wt *= bit ? w[d] : (1.0 - w[d]);
                idx += (i[d] + (bit ? 1 : 0)) * (d == 3 ? 1 : m_memStrides[d]);
            }
            res += wt * m_values[idx];
        }
		
        return res;
    }

private:
	const char* m_name;
	std::vector<std::vector<T>> m_axes;
    std::vector<T> m_values;
    std::vector<std::size_t> m_gridSizes;
    std::vector<std::size_t> m_memStrides;
    std::vector<std::size_t> m_searchStrides;
	std::vector<T> m_minValues, m_maxValues;
    int m_dim;
	
	void initializeInternal(const std::vector<std::vector<T>> &axes, const std::vector<T> &f) {
        m_axes = axes;
        m_values = f;
        m_dim = static_cast<int>(axes.size());
        precomputeStrides();
    }

	void precomputeStrides() {
		m_gridSizes.resize(m_dim);
        m_searchStrides.resize(m_dim);
        m_memStrides.assign(m_dim, 1);
		m_minValues.resize(m_dim);
        m_maxValues.resize(m_dim);

        for (int d = 0; d < m_dim; ++d) {
			const auto& vec = m_axes[d];

            std::size_t step = 1;
            while (step < m_axes[d].size() && m_axes[d][step] == m_axes[d][0]) step++;
            m_searchStrides[d] = step;

			if (vec.size() == m_values.size()) {
				std::size_t period = vec.size();
				for (std::size_t i = step; i < vec.size(); ++i) {
					if (vec[i] == vec[0] && vec[i-1] != vec[0]) {
						period = i;
						break;
					}
				}
				m_gridSizes[d] = period / step;
			} else {
				m_searchStrides[d] = 1;
                m_gridSizes[d] = vec.size();
			}

			m_minValues[d] = vec[0];
			m_maxValues[d] = vec[(m_gridSizes[d] - 1) * m_searchStrides[d]];
        }

        for (int d = m_dim - 2; d >= 0; --d) {
            m_memStrides[d] = m_gridSizes[d + 1] * m_memStrides[d + 1];
        }
	}

	inline std::size_t locateGridIndex(int axis, T val) const {
		const T minV = m_minValues[axis];
    	const T maxV = m_maxValues[axis];

		VALIDATE_BOUNDS(axis, val, minV, maxV);

		if (val <= minV) return 0;
        if (val >= maxV) return m_gridSizes[axis] - 2;

        const auto& vec = m_axes[axis];
        std::size_t ss = m_searchStrides[axis];
        std::size_t low = 0, high = m_gridSizes[axis] - 2, ans = 0;

        while (low <= high) {
            std::size_t mid = low + (high - low) / 2;
            if (vec[mid * ss] <= val) {
                ans = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return ans;
    }
};

#endif //LINEAR_INTERPOLATOR_HPP