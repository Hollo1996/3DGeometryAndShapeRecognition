#include "BezierTri.h"
#include <QGLViewer\vec.h>

using qglviewer::Vec;

//indexes of partial derivatives of Bernsteins
const int P = 0;
const int U = 1;
const int V = 2;
const int UU = 3;
const int UV = 4;
const int VV = 5;

//Constructor and Destructor
BezierTri::BezierTri(size_t degree) : degree(degree) {
	for (size_t row = 0; row <= degree; row++) {
		controls.push_back(std::vector<Vec>());
		for (size_t column = 0; column <= row; column++) {
			controls[row].push_back(Vec(0, 0, 0));
		}
	}
}

BezierTri::BezierTri() : degree(0) {}
BezierTri::~BezierTri() {}

//Reset Bezier and set new degree
void BezierTri::resize(size_t degree) {
	this->degree = degree;
	controls.clear();
	for (size_t row = 0; row <= degree; row++) {
		controls.push_back(std::vector<Vec>());
		for (size_t column = 0; column <= row; column++) {
			controls[row].push_back(Vec(0, 0, 0));
		}
	}
}
//Get degree
size_t BezierTri::getDegree() const { return degree; }

//Set all controlpoints
void BezierTri::setControlPoints(std::vector<std::vector<Vec>> newControls) {
	for (size_t row = 0; row <= degree; row++) {
		for (size_t column = 0; column <= row; column++) {
			controls[row][column] = newControls[row][column];
		}
	}
}

//Handle controlpoints by (k,l,m) discrate coordinates
void BezierTri::setControlPoint(int k, int l, Vec position) {
	controls[k + l][l] = position;
}
Vec BezierTri::getControlPoint(int k, int l) const {
	return controls[k + l][l];
}
//Handle controlpoints by (k,l,m) discrate coordinates by a (-1,1,0) directional majority index
void BezierTri::setControlPoint(int index, Vec position) {
	int row = 0;
	int column = 0;
	for (; index > row; index -= row) {
		row++;
	}
	column = index;
	controls[row][column] = position;
}
Vec BezierTri::getControlPoint(int index) const {
	int row = 0;
	int column = 0;
	for (; index > row; index -= row) {
		row++;
	}
	column = index;
	return controls[row][column];
}

//Elevate degree of Bezier
void BezierTri::elevateDegree() {
	controls.push_back(std::vector<Vec>());
	for (int column = 0; column <= degree + 1; column++) {
		controls[degree + 1].push_back(Vec(0, 0, 0));
	}

	for (int row = degree; row >= 0; row--) {
		for (int column = row; column >= 0; column--) {
			int k = row - column;
			int l = column;
			controls[row + 1][column] += controls[row][column] * (k + 1);
			controls[row + 1][column + 1] += controls[row][column] * (l + 1);
			controls[row][column] *= (degree + 1 - k - l);
		}
	}

	for (int row = 0; row <= degree + 1; row++) {
		for (int column = 0; column <= row; column++) {
			controls[row][column] /= (degree + 1);
		}
	}

	degree++;
}

//Recursive for Bernsteins and their derivatives
void BezierTri::bernsteinRecursion(double timeU, double timeV, std::vector< std::vector< std::vector< double > > >& result, int depth) {
	//Actual degree
	int actualDegree = degree - depth;
	//Setting size of result container
	if (depth == 0) {
		result.clear();

		for (int row = 0; row <= degree; row++) {
			result.push_back(std::vector< std::vector< double > >());
			for (int column = 0; column <= row; column++) {
				result[row].push_back(std::vector< double >());
				for (int der = 0; der < 6; der++) {
					result[row][column].push_back(0.0);
				}
			}
		}
		//Degree 0
		result[0][0][0] = 1;
	}
	//We allready have degree 0
	if (depth != degree) {

		//Lower degree Bernsteins
		bernsteinRecursion(timeU, timeV, result, depth + 1);

		//Second partial derivatives
		if (depth == 1) {
			for (int row = 0; row < actualDegree; row++) {
				for (int column = 0; column <= row; column++) {
					// B(U,V,W) / dU * dU
					result[row + 2][column][UU] += actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row + 1][column][UU] -= 2 * actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row][column][UU] += actualDegree * (actualDegree + 1) * result[row][column][P];

					// B(U,V,W) / dV * dV
					result[row + 2][column + 2][VV] += actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row + 1][column + 1][VV] -= 2 * actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row][column][VV] += actualDegree * (actualDegree + 1) * result[row][column][P];

					// B(U,V,W) / dU * dV
					result[row + 2][column + 1][UV] += actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row + 1][column][UV] -= actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row + 1][column + 1][UV] -= actualDegree * (actualDegree + 1) * result[row][column][P];
					result[row][column][UV] += actualDegree * (actualDegree + 1) * result[row][column][P];
				}
			}
		}

		//First partial derivatives
		if (depth == 0) {
			for (int row = 0; row < actualDegree; row++) {
				for (int column = 0; column <= row; column++) {
					// B(U,V,W) / dU
					result[row + 1][column][U] += actualDegree * result[row][column][P];
					result[row][column][U] -= actualDegree * result[row][column][P];

					// B(U,V,W) / dV
					result[row + 1][column + 1][V] += actualDegree * result[row][column][P];
					result[row][column][V] -= actualDegree * result[row][column][P];
				}
			}
		}

		//Calculating actual degree Bernsteins
		for (int row = actualDegree - 1; row >= 0; row--) {
			for (int column = row; column >= 0; column--) {
				result[row + 1][column][P] += result[row][column][P] * timeU;
				result[row + 1][column + 1][P] += result[row][column][P] * timeV;
				result[row][column][P] *= (1 - timeU - timeV);
			}
		}

	}
}

//Get all Bernstein weights and their derivatives
void BezierTri::generateBernsteinDerivatives(double timeU, double timeV, std::vector< std::vector< std::vector< double > > >& result) {
	bernsteinRecursion(timeU, timeV, result);
}
//Calculate one point from Bernstein weights and derivatives
Vec BezierTri::position(std::vector < std::vector < std::vector<double> >> derivatives) {
	Vec pos = Vec(0, 0, 0);

	for (size_t k = 0; k <= degree; k++)
		for (size_t l = 0; l <= k; l++)
			pos += controls[k][l] * derivatives[k][l][P];

	return pos;
}

//Calculate normalvector in one position from Bernstein weights and derivatives
Vec BezierTri::normal(std::vector < std::vector < std::vector<double> > > derivatives) {
	Vec du(0.0, 0.0, 0.0);
	Vec dv(0.0, 0.0, 0.0);

	for (size_t k = 0; k <= degree; k++)
		for (size_t l = 0; l <= k; l++) {
			du += controls[k][l] * derivatives[k][l][U];
			dv += controls[k][l] * derivatives[k][l][V];
		}

	Vec cross = (du ^ dv);
	Vec normal = cross / sqrt(cross * cross);
	return normal;
}

//Calculate mean Curvature in one position from Bernstein weights and derivatives
float BezierTri::meanCurvature(std::vector < std::vector < std::vector<double> > > derivatives) {
	float mean;
	Vec du(0.0, 0.0, 0.0);
	Vec dv(0.0, 0.0, 0.0);
	Vec duu(0.0, 0.0, 0.0);
	Vec dvv(0.0, 0.0, 0.0);
	Vec duv(0.0, 0.0, 0.0);

	for (size_t k = 0; k <= degree; k++)
		for (size_t l = 0; l <= k; l++) {
			du += controls[k][l] * derivatives[k][l][U];
			dv += controls[k][l] * derivatives[k][l][V];
			duu += controls[k][l] * derivatives[k][l][UU];
			dvv += controls[k][l] * derivatives[k][l][VV];
			duv += controls[k][l] * derivatives[k][l][UV];
		}

	double E = du * du;
	double F = du * dv;
	double G = dv * dv;
	Vec normalDir = normal(derivatives);
	double L = normalDir * duu;
	double M = normalDir * duv;
	double N = normalDir * dvv;
	mean = (N * E - 2 * M * F + L * G) / (2 * (E * G - F * F));

	return mean;
}

Vec BezierTri::updateAlfaBeta(Vec common0, Vec common1, Vec master, Vec slave) {
	Vec c0m = master - common0;
	Vec c0s = slave - common0;
	Vec ms = slave - master;
	Vec c1m = master - common1;
	Vec c1s = slave - common1;

	double c0mLength = sqrt(c0m * c0m);
	double c0sLength = sqrt(c0s * c0s);
	double msLength = sqrt(ms * ms);
	double c1mLength = sqrt(c1m * c1m);
	double c1sLength = sqrt(c1s * c1s);

	double s0 = (c0mLength + c0sLength + msLength) / 2;
	double s1 = (c1mLength + c1sLength + msLength) / 2;

	double T0 = sqrt(s0 * (s0 - c0mLength) * (s0 - c0sLength) * (s0 - msLength));
	double T1 = sqrt(s1 * (s1 - c1mLength) * (s1 - c1sLength) * (s1 - msLength));

	double beta = T0 / (T0 + T1);
	Vec intersection = common0 + (common1 - common0) * beta;
	double alfa = sqrt((intersection - master) * (intersection - master)) / msLength;
	return Vec(alfa, beta, 0);
}

void BezierTri::updateSlave(BezierTri& master, double alfa0, double alfa1, double beta0, double beta1) {


	std::vector<Vec> masterParalell = std::vector<Vec>();
	std::vector<Vec> masterCommon = std::vector<Vec>();
	std::vector<Vec> slaveCommon = std::vector<Vec>();
	std::vector<Vec> slaveParalell = std::vector<Vec>();

	for (int column = 0; column < degree; column++) {
		masterCommon.push_back(master.controls[degree][degree - column]);
		masterParalell.push_back(master.controls[degree - 1][degree - 1 - column]);
		slaveCommon.push_back(masterCommon[column]);
		slaveParalell.push_back(controls[degree - 1][column]);
	}

	masterCommon.push_back(master.controls[degree][0]);
	slaveCommon.push_back(masterCommon[degree]);


	slaveParalell[0] = (
		beta0 * masterCommon[0]
		+ (1 - beta0) * masterCommon[1]
		- alfa0 * masterParalell[0]
		) / (1 - alfa0);
	slaveParalell[degree - 1] = (
		beta1 * masterCommon[degree - 1]
		+ (1 - beta1) * masterCommon[degree]
		- alfa1 * masterParalell[degree - 1]
		) / (1 - alfa1);

	for (int column = 1; column < degree-1; column++) {
		slaveParalell[column] = (
			((double)column / degree) *
			( beta1 * masterCommon[column - 1]
				+ (1 - beta1) * masterCommon[column]
				- alfa1 * masterParalell[column - 1]
				- (1 - alfa1) * slaveParalell[column - 1]
				)
			 + (1 - (double)column / degree) * (
				beta0 * masterCommon[column]
				+ (1 - beta0) * masterCommon[column + 1]
				- alfa0 * masterParalell[column]
				)
			) / ((1 - (double)column / degree) * (1 - alfa0));
	}


	for (int column = 0; column < degree; column++) {
		controls[degree][column] = slaveCommon[column];
		controls[degree - 1][column] = slaveParalell[column];
	}
	controls[degree][degree] = slaveCommon[degree];
}
