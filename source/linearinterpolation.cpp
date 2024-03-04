#include "linearinterpolation.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

//CONSTRUCTORS:
interpolationF::interpolationF() {

}

//input is 2 1D arrays:
interpolationF::interpolationF(const double *xData, const double *fData, size_t NofElements)
{
	SetData(xData, fData, NofElements);
}

void interpolationF::SetData(const double *xData, const double *fData, size_t NofElements)
{
	m_variableN = 1;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(xData, xData + m_dataLength);
    m_data[1] = std::vector<double>(fData, fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 2 1D vectors:
interpolationF::interpolationF(const std::vector<double> &xData, const std::vector<double> &fData)
{
	SetData(xData, fData);
}

void interpolationF::SetData(const std::vector<double> &xData, const std::vector<double> &fData)
{
	m_variableN = 1;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(xData.begin(), xData.begin() + m_dataLength);
    m_data[1] = std::vector<double>(fData.begin(), fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 3 1D arrays:
interpolationF::interpolationF(const double *x1Data, const double *x2Data, const double *fData, size_t NofElements)
{
	SetData(x1Data, x2Data, fData, NofElements);
}

void interpolationF::SetData(const double *x1Data, const double *x2Data, const double *fData, size_t NofElements)
{
	m_variableN = 2;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<double>(x2Data, x2Data + m_dataLength);
    m_data[2] = std::vector<double>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 3 1D vectors:
interpolationF::interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &fData)
{
	SetData(x1Data, x2Data, fData);
}

void interpolationF::SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &fData)
{
	m_variableN = 2;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<double>(x2Data.begin(), x2Data.begin() + m_dataLength);
    m_data[2] = std::vector<double>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 2 1D vectors (grids) and 1 2d vector (function values):
interpolationF::interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<std::vector<double>> &fData)
{
	SetData(x1Data, x2Data, fData);
}

void interpolationF::SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<std::vector<double>> &fData)
{
	m_variableN = 2;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data.begin(), x1Data.end());
	m_data[1] = std::vector<double>(x2Data.begin(), x2Data.end());
	for (const auto &row : fData)
		for (const auto &elem : row)
			m_data[2].push_back(elem);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 4 1D arrays:
interpolationF::interpolationF(const double *x1Data, const double *x2Data, const double *x3Data, const double *fData, size_t NofElements)
{
	SetData(x1Data, x2Data, x3Data, fData, NofElements);
}

void interpolationF::SetData(const double *x1Data, const double *x2Data, const double *x3Data, const double *fData, size_t NofElements)
{
	m_variableN = 3;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<double>(x2Data, x2Data + m_dataLength);
	m_data[2] = std::vector<double>(x3Data, x3Data + m_dataLength);
    m_data[3] = std::vector<double>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 4 1D vectors:
interpolationF::interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &fData)
{
	SetData(x1Data, x2Data, x3Data, fData);
}

void interpolationF::SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &fData)
{
	m_variableN = 3;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<double>(x2Data.begin(), x2Data.begin() + m_dataLength);
	m_data[2] = std::vector<double>(x3Data.begin(), x3Data.begin() + m_dataLength);
    m_data[3] = std::vector<double>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 5 1D arrays:
interpolationF::interpolationF(const double *x1Data, const double *x2Data, const double *x3Data, const double *x4Data, const double *fData, size_t NofElements)
{
	SetData(x1Data, x2Data, x3Data, x4Data, fData, NofElements);
}

void interpolationF::SetData(const double *x1Data, const double *x2Data, const double *x3Data, const double *x4Data, const double *fData, size_t NofElements)
{
	m_variableN = 4;
	m_dataLength = NofElements;

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data, x1Data + m_dataLength);
	m_data[1] = std::vector<double>(x2Data, x2Data + m_dataLength);
	m_data[2] = std::vector<double>(x3Data, x3Data + m_dataLength);
	m_data[3] = std::vector<double>(x4Data, x4Data + m_dataLength);
    m_data[4] = std::vector<double>( fData,  fData + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//input is 5 1D vectors:
interpolationF::interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &x4Data, const std::vector<double> &fData)
{
	SetData(x1Data, x2Data, x3Data, x4Data, fData);
}

void interpolationF::SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &x4Data, const std::vector<double> &fData)
{
	m_variableN = 4;
	m_dataLength = fData.size();

	m_data.resize(m_variableN+1);

	m_data[0] = std::vector<double>(x1Data.begin(), x1Data.begin() + m_dataLength);
	m_data[1] = std::vector<double>(x2Data.begin(), x2Data.begin() + m_dataLength);
	m_data[2] = std::vector<double>(x3Data.begin(), x3Data.begin() + m_dataLength);
	m_data[3] = std::vector<double>(x4Data.begin(), x4Data.begin() + m_dataLength);
    m_data[4] = std::vector<double>( fData.begin(),  fData.begin() + m_dataLength);

	createGrids();

	for (size_t iv=0; iv<m_variableN; iv++)
		if (m_data[iv].size() < 2)
			std::cerr << "Error: not enough data for interplation for variable " + std::to_string(iv) + "." << std::endl;
}

//DESTRUCTORS:
interpolationF::~interpolationF() {

}

//INTERPOLATION FUNCTIONS:
//1D interpolation
double interpolationF::interpolation(double pointValue) const
{
	if (m_variableN > 1) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else {
		if (pointValue < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		return interpolation1D(pointValue);
	}
}

//2D interpolation
double interpolationF::interpolation(double pointValue1, double pointValue2) const
{
	if (m_variableN < 2) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if (m_variableN > 2) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		return interpolation2D(pointValue1, pointValue2);
	}
}

//3D interpolation
double interpolationF::interpolation(double pointValue1, double pointValue2, double pointValue3) const
{
	if (m_variableN < 3) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if (m_variableN > 3) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue3 < m_domain[2][0]) {
			std::cerr << "Error: point value in dimension 3 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue3 > m_domain[2][1]) {
			std::cerr << "Error: point value in dimension 3 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		return interpolation3D(pointValue1, pointValue2, pointValue3);
	}
	return 0.0;
}

//4D interpolation
double interpolationF::interpolation(double pointValue1, double pointValue2, double pointValue3, double pointValue4) const
{
	if (m_variableN < 4) {
		std::cerr << "Error: too much points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else if (m_variableN > 4) {
		std::cerr << "Error: not enough points for interpolation." << std::endl;
		return std::numeric_limits<double>::quiet_NaN();
	}
	else {
		if (pointValue1 < m_domain[0][0]) {
			std::cerr << "Error: point value in dimension 1 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue1 > m_domain[0][1]) {
			std::cerr << "Error: point value in dimension 1 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 < m_domain[1][0]) {
			std::cerr << "Error: point value in dimension 2 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue2 > m_domain[1][1]) {
			std::cerr << "Error: point value in dimension 2 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue3 < m_domain[2][0]) {
			std::cerr << "Error: point value in dimension 3 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue3 > m_domain[2][1]) {
			std::cerr << "Error: point value in dimension 3 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue4 < m_domain[3][0]) {
			std::cerr << "Error: point value in dimension 4 smaller than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		if (pointValue4 > m_domain[3][1]) {
			std::cerr << "Error: point value in dimension 4 larger than domain." << std::endl;
			return std::numeric_limits<double>::quiet_NaN();
		}
		return interpolation4D(pointValue1, pointValue2, pointValue3, pointValue4);
	}
	return 0.0;
}

const std::vector<std::vector<double>> & interpolationF::domain() const
{
	return m_domain;
}

const std::vector<double> & interpolationF::codomain() const
{
	return m_codomain;
}

void interpolationF::createGrids()
{
	for (size_t iv=0; iv<m_variableN; iv++) {
		std::sort(m_data[iv].begin(), m_data[iv].end());
		m_data[iv].erase(std::unique(m_data[iv].begin(), m_data[iv].end()), m_data[iv].end());
		m_domain.push_back({m_data[iv].front(), m_data[iv].back()});
	}
	m_codomain.push_back(*std::min_element(m_data[m_variableN].begin(), m_data[m_variableN].end()));
	m_codomain.push_back(*std::max_element(m_data[m_variableN].begin(), m_data[m_variableN].end()));
}

// void interpolationF::locatePointF(const std::vector<double> &points, std::vector<size_t> &positions) const
// {
// 	positions.resize(points.size(), 0.0L);
// 	size_t left, right, middle;
// 	for (size_t ip=0; ip<points.size(); ip++) {
// 		left = 0;
// 		right = m_data[ip].size() - 1;
// 		positions[ip] = left;
// 		while (left <= right) {
// 			middle = left + (right - left) / 2;
// 			if (m_data[ip][middle] < points[ip]) {
// 				left = middle + 1;
// 			}
// 			else if (m_data[ip][middle] > points[ip]) {
// 				right = middle - 1;
// 			}
// 			else {
// 				positions[ip] = middle;
// 				break;
// 			}
// 			if (std::abs(m_data[ip][middle] - points[ip]) < std::abs(m_data[ip][positions[ip]] - points[ip])) {
// 				positions[ip] = middle;
// 			}
// 		}
// 	}
// }

void interpolationF::locatePointF(const std::vector<double> &points, std::vector<size_t> &positions) const
{
	positions.resize(points.size(), 0);
    int ju, jm, jl, mm = 1 + 1;
	bool ascnd;
	for (size_t iv=0; iv<m_data.size()-1; iv++)
	{
		jl = 0;
		ju = m_data[iv].size() - 1;
		ascnd = (m_data[iv].back() >= m_data[iv][0]);

		while ((ju - jl) >1) {
			jm = (ju + jl) >> 1;
			if ((points[iv] >= m_data[iv][jm]) == ascnd) {
				jl = jm;
			}
			else {
				ju = jm;
			}
		}
		int n = static_cast<int>(m_data[iv].size());
		positions[iv] = static_cast<size_t>(std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1))));
	}
}

//1D linear interpolation
double interpolationF::lin1DInterpolation(const double x[2], const double f[2], double xx) const
{
	return (f[0] + (xx - x[0])*(f[1] - f[0]) / (x[1] - x[0]));
}

//1D interpolation (full function)
double interpolationF::interpolation1D(double pointValue) const
{
	//searching for position
	const std::vector<double> points{pointValue};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	//setting x and Q values
	double x[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	double Q[] = {m_data[1][positions[0]], m_data[1][positions[0] + 1]};

	return lin1DInterpolation(x, Q, pointValue);
}

//2D interpolation
double interpolationF::interpolation2D(double pointValue1, double pointValue2) const
{
	//searching for position
	const std::vector<double> points{pointValue1, pointValue2};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	double x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	double x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};	

	double Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = m_data[2][(positions[0] + i1)*m_data[1].size() + (positions[1] + i2)];

	double Q1[2];
	for (int i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}

//3D interpolation
double interpolationF::interpolation3D(double pointValue1, double pointValue2, double pointValue3) const
{
	//searching for position
	const std::vector<double> points{pointValue1, pointValue2, pointValue3};
	std::vector<size_t> positions;
	locatePointF(points, positions);
	
	double x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	double x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};
	double x3[] = {m_data[2][positions[2]], m_data[2][positions[2] + 1]};

	double Q3[2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				Q3[i1][i2][i3] = m_data[3][(positions[0] + i1)*m_data[2].size()*m_data[1].size() + 
										   (positions[1] + i2)*m_data[2].size() +
										   (positions[2] + i3)]; 

	double Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = lin1DInterpolation(x3, Q3[i1][i2], pointValue3);

	double Q1[2];
	for (size_t i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}

//4D interpolation
double interpolationF::interpolation4D(double pointValue1, double pointValue2, double pointValue3, double pointValue4) const
{
	//searching for position
	const std::vector<double> points{pointValue1, pointValue2, pointValue3, pointValue4};
	std::vector<size_t> positions;
	locatePointF(points, positions);

	double x1[] = {m_data[0][positions[0]], m_data[0][positions[0] + 1]};
	double x2[] = {m_data[1][positions[1]], m_data[1][positions[1] + 1]};
	double x3[] = {m_data[2][positions[2]], m_data[2][positions[2] + 1]};
	double x4[] = {m_data[3][positions[3]], m_data[3][positions[3] + 1]};

	double Q4[2][2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				for (size_t i4=0; i4<2; i4++)
					Q4[i1][i2][i3][i4] = m_data[4][(positions[0] + i1)*m_data[3].size()*m_data[2].size()*m_data[1].size() +
												   (positions[1] + i2)*m_data[3].size()*m_data[2].size() +
												   (positions[2] + i3)*m_data[3].size() +
												   (positions[3] + i4)];

	double Q3[2][2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			for (size_t i3=0; i3<2; i3++)
				Q3[i1][i2][i3] = lin1DInterpolation(x4, Q4[i1][i2][i3], pointValue4);
	
	double Q2[2][2];
	for (size_t i1=0; i1<2; i1++)
		for (size_t i2=0; i2<2; i2++)
			Q2[i1][i2] = lin1DInterpolation(x3, Q3[i1][i2], pointValue3);
	
	
	double Q1[2];
	for (size_t i1=0; i1<2; i1++)
		Q1[i1] = lin1DInterpolation(x2, Q2[i1], pointValue2);

	return lin1DInterpolation(x1, Q1, pointValue1);
}