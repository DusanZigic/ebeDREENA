#ifndef HEADERFILE_LINEARINTERPOLATION
#define HEADERFILE_LINEARINTERPOLATION

#include <vector>

class interpolationF {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//public functions:
public:

	//CONSTRUCTORS:

	interpolationF();

	//input is 2 1D arrays:
	interpolationF(const double *xData, const double *fData, size_t NofElements);
	void SetData(const double *xData, const double *fData, size_t NofElements);

	//input is 2 1D vectors:
	interpolationF(const std::vector<double> &xData, const std::vector<double> &fData);
	void SetData(const std::vector<double> &xData, const std::vector<double> &fData);

	//input is 3 1D arrays:
	interpolationF(const double *x1Data, const double *x2Data, const double *fData, size_t NofElements);
	void SetData(const double *x1Data, const double *x2Data, const double *fData, size_t NofElements);

	//input is 3 1D vectors:
	interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &fData);
	void SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &fData);

	//input is 2 1D vectors (grids) and 1 2d vector (function values):
	interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<std::vector<double>> &fData);
	void SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<std::vector<double>> &fData);

	//input is 4 1D arrays:
	interpolationF(const double *x1Data, const double *x2Data, const double *x3Data, const double *fData, size_t NofElements);
	void SetData(const double *x1Data, const double *x2Data, const double *x3Data, const double *fData, size_t NofElements);

	//input is 4 1D vectors:
	interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &fData);
	void SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &fData);

	//input is 5 1D arrays:
	interpolationF(const double *x1Data, const double *x2Data, const double *x3Data, const double *x4Data, const double *fData, size_t NofElements);
	void SetData(const double *x1Data, const double *x2Data, const double *x3Data, const double *x4Data, const double *fData, size_t NofElements);

	//input is 5 1D vectors:
	interpolationF(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &x4Data, const std::vector<double> &fData);
	void SetData(const std::vector<double> &x1Data, const std::vector<double> &x2Data, const std::vector<double> &x3Data, const std::vector<double> &x4Data, const std::vector<double> &fData);

	//DESTRUCTOR:
	~interpolationF();

	//INTERPOLATION FUNCTIONS:
	//1D interpolation
	double interpolation(double pointValue) const;

	//2D interpolation
	double interpolation(double pointValue1, double pointValue2) const;

	//3D interpolation
	double interpolation(double pointValue1, double pointValue2, double pointValue3) const;

	//4D interpolation
	double interpolation(double pointValue1, double pointValue2, double pointValue3, double pointValue4) const;

	//miscellaneous FUNCTIONS:

	//function that returns domains:
	const std::vector<std::vector<double>> & domain() const;

	//function that returns codomain:
	const std::vector<double> & codomain() const;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//private variables and functions:
private:

	size_t m_dataLength;
	std::vector<std::vector<double>> m_data;
	size_t m_variableN;
	std::vector<size_t> m_gridLengths;
	std::vector<size_t> m_relPosition;
	std::vector<std::vector<double>> m_domain;
	std::vector<double> m_codomain;

	void createGrids();

	//function that locates points
	void locatePointF(const std::vector<double> &points, std::vector<size_t> &positions) const;

	double lin1DInterpolation(const double x[2], const double f[2], double xx) const;

	//1D interpolation (full function)
	double interpolation1D(double pointValue) const;

	//2D interpolation
	double interpolation2D(double pt1, double pt2) const;

	//3D interpolation
	double interpolation3D(double pt1, double pt2, double pt3) const;

	//4D interpolation
	double interpolation4D(double pt1, double pt2, double pt3, double pt4) const;
};

#endif