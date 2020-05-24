#pragma once
#include <string>
#include <QGLViewer\vec.h>
#include <BezierRect.h>

using qglviewer::Vec;


//Class for storing and handleing Bezier Triangles
class BezierTri : public QObject
{
	Q_OBJECT
private:
	//degree of Bernstein
	size_t degree;
	//Control points. First row has 1 value, second has 2 ...
	std::vector< std::vector< Vec > > controls;
	//Recursive for Bernsteins and their derivatives
	void bernsteinRecursion(double timeU, double timeV, std::vector< std::vector< std::vector< double > > >& result, int depth = 0);
public:
	//Constructor and Destructor
	BezierTri();
	BezierTri(size_t degree);
	~BezierTri();

	//Reset Bezier and set new degree
	void resize(size_t degree);
	//Get degree
	size_t getDegree() const;

	//Set all controlpoints
	void setControlPoints(std::vector<std::vector<Vec>> newControls);
	//Handle controlpoints by (k,l,m) discrate coordinates
	void setControlPoint(int k, int l, Vec position);
	Vec getControlPoint(int k, int l) const;
	//Handle controlpoints by (k,l,m) discrate coordinates by a (-1,1,0) directional majority index
	void setControlPoint(int index, Vec position);
	Vec getControlPoint(int index) const;

	//Get all Bernstein weights and their derivatives
	void generateBernsteinDerivatives(double timeU, double timeV, std::vector < std::vector<std::vector<double>> >& result);
	//Calculate one point from Bernstein weights and derivatives
	Vec position(std::vector < std::vector < std::vector<double> >> derivatives);
	//Calculate normalvector in one position from Bernstein weights and derivatives
	Vec normal(std::vector < std::vector < std::vector<double> > >derivatives);
	//Calculate mean Curvature in one position from Bernstein weights and derivatives
	float meanCurvature(std::vector < std::vector < std::vector<double> >> derivatives);

	//Elevate degree of Bezier
	void elevateDegree();

	Vec updateAlfaBeta(Vec common0, Vec common1, Vec master, Vec slave);

	void updateSlave(BezierTri& master, double alfa0, double alfa1, double beta0, double beta1);


	//Not Implemented!!!
	//Elevate degree of a Bezier curve
	//void elevateCurveDegree(std::vector < Vec > controlPonóints);
	//Transform Bezier triangle to Bezier Rectangle
	//BezierRect transform();
};

