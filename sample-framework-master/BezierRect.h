#pragma once
#include <vector>
#include <string>
#include <QGLViewer\vec.h>

using qglviewer::Vec;

//Class for storing and handleing Bezier Rectangels
class BezierRect : public QObject
{
	Q_OBJECT
	private:
		//Degrees of Bernsteins
		size_t degreeK, degreeL;
		//Controlpoints
		std::vector<std::vector<Vec>> controls;
		//Recursive function for Bernstein and its primary and secundary derivatives.
		void bernsteinRecursion(double time, std::vector<Vec>& result, size_t maxLevel, int depth = 0);
	public:
		// constructor and destructor
		BezierRect(size_t degreeK, size_t degreeL);
		~BezierRect();

		//Resets bezier and sets given Bernstein levels
		void resize(size_t degreeK, size_t degreeL);
		//Get levels of Bernsteins
		size_t getK() const;
		size_t getL() const;

		//Set all control points
		void setControlPoints(std::vector<std::vector<Vec>> newControls);
		//Set one conntrol point given by discrate coordinates.
		void setControlPoint(int k, int l, Vec position);
		Vec getControlPoint(int k, int l) const;
		//Set one conntrol point given by k major index.
		void setControlPoint(int kMajorIndex, Vec position);
		Vec getControlPoint(int kMajorIndex) const;

		//Generate all Bernsteins and their derivatives
		void generateBernsteinDerivatives( double timeU, double timeV, std::vector<Vec>& resultU, std::vector<Vec>& resultV);
		//Calculate one point from Bernstein weights and derivatives
		Vec position(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV);
		//Calculate normalvector in one position from Bernstein weights and derivatives
		Vec normal(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV);
		//Calculate mean Curvature in one position from Bernstein weights and derivatives
		float meanCurvature(std::vector < Vec > derivativesU, std::vector < Vec > derivativesV);
		
		//Elevate degree of Bezier
		void elevateDegree();
};

