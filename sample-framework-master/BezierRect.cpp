#include "BezierRect.h"
#include <QGLViewer\vec.h>

using qglviewer::Vec;

//Constructor
BezierRect::BezierRect(size_t degreeK, size_t degreeL) : degreeK(degreeK), degreeL(degreeL) {
	for (size_t i = 0; i <= degreeK; i++) {
		controls.push_back(std::vector<Vec>());
		for (size_t j = 0; j <= degreeL; j++) {
			controls[i].push_back(Vec(0, 0, 0));
		}
	}
}

//Destructor
BezierRect::~BezierRect() {}

//Reset Bezier and set new degrees
void BezierRect::resize(size_t degreeK, size_t degreeL) {
	this->degreeK = degreeK;
	this->degreeL = degreeL;
	controls.clear();
	for (size_t k = 0; k <= degreeK; k++) {
		controls.push_back(std::vector<Vec>());
		for (size_t l = 0; l <= degreeL; l++) {
			controls[k].push_back(Vec(0, 0, 0));
		}
	}
}

//Get degrees
size_t BezierRect::getK() const { return degreeK; }

size_t BezierRect::getL() const { return degreeL; }

//Set all controlpoints
void BezierRect::setControlPoints(std::vector<std::vector<Vec>> newControls) {
	for (size_t k = 0; k <= degreeK; k++) {
		for (size_t l = 0; l <= degreeL; l++) {
			controls[k][l] = newControls[k][l];
		}
	}
}

//set and get control points by descreet coordinates
void BezierRect::setControlPoint(int k, int l, Vec position) {
	controls[k][l] = position;
}
Vec BezierRect::getControlPoint(int k, int l) const {
	return controls[k][l];
}
//set and get control points by k major index
void BezierRect::setControlPoint(int rowMajorIndex, Vec position) {
	controls[rowMajorIndex/(degreeL+1)][rowMajorIndex % (degreeL+1)] = position;
}
Vec BezierRect::getControlPoint(int rowMajorIndex) const {
	return controls[rowMajorIndex / (degreeL + 1)][rowMajorIndex % (degreeL + 1)];
}

//elevate degree of Bezier
void BezierRect::elevateDegree() {

	//Eleveta K degree
	std::vector<std::vector<Vec>> temporalControlPoints;
	float rowCount = degreeK + 1;
	float columnCount = degreeL + 1;
	for (size_t row = 0; row < rowCount; row++) {
		temporalControlPoints.push_back(std::vector<Vec>());
		temporalControlPoints[row].push_back(controls[row][0]);
		for (size_t column = 1; column < columnCount; column++) {
			temporalControlPoints[row].push_back(
				(controls[row][column - 1] * column / columnCount)
				+ (controls[row][column] * (columnCount - column) / columnCount)
			);
		}
		temporalControlPoints[row].push_back(controls[row][columnCount - 1]);
	}

	//set new size
	resize(++degreeK, ++degreeL);

	//Eleveta L degree
	//Copy first row
	for (size_t column = 0; column < columnCount + 1; column++) {
		controls[0][column] = temporalControlPoints[0][column];
	}
	//Calculate intermadiet rows
	for (size_t row = 1; row < rowCount; row++) {
		for (size_t column = 0; column < columnCount + 1; column++) {
			controls[row][column] =
				temporalControlPoints[(row - 1)][column] * row / rowCount
				+ temporalControlPoints[row][column] * (rowCount - row) / rowCount;
		}
	}
	//Copy last row
	for (size_t column = 0; column < columnCount + 1; column++) {
		controls[rowCount][column] = temporalControlPoints[rowCount-1][column];
	}
}

//Calculation of all Bernstein weights and their derivatives
void BezierRect::bernsteinRecursion(double time, std::vector<Vec>& result, size_t maxdegree, int depth) {

	//Actual degree
	int degree = maxdegree - depth;

	//In first round of recursion we set the result vector's size
	if (depth == 0) {
		result.clear();
		result.reserve(degree + 1);
		for (int i = 0; i <= degree; i++) {
			result.push_back(Vec(0.0, 0.0, 0.0));
		}
		//Degree 0
		result[0] = Vec(1.0, 0.0, 0.0);
	}
	//We already have the result of degree 0
	if (degree != 0) {

		double t1 = 1.0 - time;

		//Calculating lower degree
		bernsteinRecursion(time, result, maxdegree, depth + 1);
		
		//Second partial derivatives
		if (depth == 1) {
			for (int i = 0; i < degree; i++) {
				result[i][2] += degree * (degree + 1) * result[i][0];
				result[i + 1][2] -= degree * (degree + 1) * result[i][0] * 2;
				result[i + 2][2] += degree * (degree + 1) * result[i][0];
			}
		}

		//Calculate first partial derivatives
		if (depth == 0) {
			for (int i = 0; i < degree; i++) {
				result[i][1] -= degree * result[i][0];
				result[i + 1][1] += degree * result[i][0];
			}
		}

		//calculate Bernstein for actual degree
		for (int i = degree - 1; i >= 0; i--) {
			result[i + 1][0] += time * result[i][0];
			result[i][0] *= t1;
		}
	}
}

//Generate all Bernsteins and their derivatives
void BezierRect::generateBernsteinDerivatives(double timeU, double timeV, std::vector<Vec>& resultU, std::vector<Vec>& resultV) {
	bernsteinRecursion(timeU, resultU, degreeK);
	bernsteinRecursion(timeV, resultV, degreeL);
}

//Calculate one point from Bernstein weights and derivatives
Vec BezierRect::position(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV) {
	Vec pos = Vec(0, 0, 0);

	for (int k = 0; k <= degreeK; k++)
		for (int l = 0; l <= degreeL; l++)
			pos += controls[k][l] * derivativesU[k][0] * derivativesV[l][0];

	return pos;
}

//Calculate normalvector in one position from Bernstein weights and derivatives
Vec BezierRect::normal(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV) {
	//Partial derivatives
	Vec du(0.0, 0.0, 0.0);
	Vec dv(0.0, 0.0, 0.0);

	//Sum
	for (int k = 0; k <= degreeK; k++)
		for (int l = 0; l <= degreeL; l++) {
			du += controls[k][l] * derivativesU[k][1] * derivativesV[l][0];
			dv += controls[k][l] * derivativesU[k][0] * derivativesV[l][1];
		}

	//calculate normal
	Vec cross = (du ^ dv);
	Vec normal = cross / sqrt(cross * cross);
	return normal;
}

//Calculate mean Curvature in one position from Bernstein weights and derivatives
float BezierRect::meanCurvature(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV) {
	//Partial derivatives
	Vec du(0.0, 0.0, 0.0);
	Vec dv(0.0, 0.0, 0.0);
	Vec duu(0.0, 0.0, 0.0);
	Vec dvv(0.0, 0.0, 0.0);
	Vec duv(0.0, 0.0, 0.0);

	//Sum
	for (int k = 0; k <= degreeK; k++)
		for (int l = 0; l <= degreeL; l++) {
			du += controls[k][l] * derivativesU[k][1] * derivativesV[l][0];
			dv += controls[k][l] * derivativesU[k][0] * derivativesV[l][1];
			duu += controls[k][l] * derivativesU[k][2] * derivativesV[l][0];
			dvv += controls[k][l] * derivativesU[k][0] * derivativesV[l][2];
			duv += controls[k][l] * derivativesU[k][1] * derivativesV[l][1];
		}

	//Partial results
	double E = du * du;
	double F = du * dv;
	double G = dv * dv;
	Vec normalDir = normal(derivativesU, derivativesV);
	double L = normalDir * duu;
	double M = normalDir * duv;
	double N = normalDir * dvv;

	//Result
	float mean = (N * E - 2 * M * F + L * G) / (2 * (E * G - F * F));

	return mean;
}
