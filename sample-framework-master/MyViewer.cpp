#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <QtGui/QKeyEvent>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

// #define BETTER_MEAN_CURVATURE

#ifdef BETTER_MEAN_CURVATURE
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#endif

#include "MyViewer.h"

#ifdef _WIN32
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_BGRA 0x80E1
#endif

MyViewer::MyViewer(QWidget* parent) :
	QGLViewer(parent), model_type(ModelType::NONE),
	mean_min(0.0), mean_max(0.0), cutoff_ratio(0.05),
	show_control_points(true), show_solid(true), show_wireframe(false),
	visualization(Visualization::PLAIN), slicing_dir(0, 0, 1), slicing_scaling(1),
	last_filename(""),
	bezierRect(0, 0)
{
	setSelectRegionWidth(10);
	setSelectRegionHeight(10);
	axes.shown = false;
}

MyViewer::~MyViewer() {
	glDeleteTextures(1, &isophote_texture);
	glDeleteTextures(1, &environment_texture);
	glDeleteTextures(1, &slicing_texture);
}

//Loading a Bezier Rectangle from file
bool MyViewer::openBezierRect(const std::string& filename, bool update_view) {
	size_t n, m;
	try {
		std::ifstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f >> n >> m;
		bezierRect.resize(n, m);
		for (size_t k = 0; k <= bezierRect.getK(); k++)
			for (size_t l = 0; l <= bezierRect.getL(); l++) {
				Vec p;
				f >> p[0] >> p[1] >> p[2];
				bezierRect.setControlPoint(k, l, p);
			}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	model_type = ModelType::BEZIER_SURFACE_RECTANGULAR;
	last_filename = filename;
	updateMesh(update_view);
	if (update_view)
		setupCamera();
	return true;
}

//Loading a Bezier Triangle from file
bool MyViewer::openBezierTri(const std::string& filename, bool update_view) {
	size_t n;
	try {
		std::ifstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f >> n;
		bezierTri[Master].resize(n);
		for (size_t row = 0; row <= bezierTri[Master].getDegree(); row++)
			for (size_t column = 0; column <= row; column++) {
				int k = row - column;
				int l = column;
				Vec p;
				f >> p[0] >> p[1] >> p[2];
				bezierTri[Master].setControlPoint(k, l, p);
			}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	model_type = ModelType::BEZIER_SURFACE_TRIANGULAR;
	last_filename = filename;
	updateMesh(update_view);
	if (update_view)
		setupCamera();
	return true;
}
bool MyViewer::openBezierTriPair(const std::string& filename, bool update_view) {
	size_t n;
	try {
		std::ifstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f >> n;
		for (int index = 0; index < 2; index++) {
			bezierTri[index].resize(n);
			for (size_t row = 0; row <= bezierTri[index].getDegree(); row++)
				for (size_t column = 0; column <= row; column++) {
					int k = row - column;
					int l = column;
					Vec p;
					f >> p[0] >> p[1] >> p[2];
					bezierTri[index].setControlPoint(k, l, p);
				}
		}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	model_type = ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR;
	last_filename = filename;
	updateMesh(update_view);
	if (update_view)
		setupCamera();
	return true;
}

//Saving a Bezier Rectangle to file
bool MyViewer::saveBezierRect(const std::string& filename) {
	if (model_type != ModelType::BEZIER_SURFACE_RECTANGULAR)
		return false;

	try {
		std::ofstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f << bezierRect.getK() << ' ' << bezierRect.getL() << std::endl;
		for (int k = 0; k <= bezierRect.getK(); k++)
			for (int l = 0; l <= bezierRect.getL(); l++) {
				Vec p = bezierRect.getControlPoint(k, l);
				f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
			}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	return true;
}

//Save a Bezier Triangle to file
bool MyViewer::saveBezierTri(const std::string& filename) {
	if (model_type != ModelType::BEZIER_SURFACE_TRIANGULAR)
		return false;

	try {
		std::ofstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f << bezierTri[Master].getDegree() << std::endl;
		for (size_t row = 0; row <= bezierTri[Master].getDegree(); row++)
			for (size_t column = 0; column <= row; column++) {
				int k = row - column;
				int l = column;
				Vec p = bezierTri[Master].getControlPoint(k, l);
				f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
			}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	return true;
}

//Save a Bezier Triangle to file
bool MyViewer::saveBezierTriPair(const std::string& filename) {
	if (model_type != ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR)
		return false;

	try {
		std::ofstream f(filename.c_str());
		f.exceptions(std::ios::failbit | std::ios::badbit);
		f << bezierTri[Master].getDegree() << std::endl;
		for (int index = 0; index < 2; index++)
			for (size_t row = 0; row <= bezierTri[index].getDegree(); row++)
				for (size_t column = 0; column <= row; column++) {
					int k = row - column;
					int l = column;
					Vec p = bezierTri[index].getControlPoint(k, l);
					f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
				}
	}
	catch (std::ifstream::failure&) {
		return false;
	}
	return true;
}

//Generate Bezier Rectangle Mesh and Mean curvature
void MyViewer::generateBezierRectMesh() {
	size_t resolution = 50;

	mesh.clear();
	std::vector<MyMesh::VertexHandle> handles, tri;
	std::vector< float > means;

	std::vector<Vec> bernsteinDerivativesU = std::vector<Vec>();
	std::vector<Vec> bernsteinDerivativesV = std::vector<Vec>();

	//Calculate vertexes and their curvature
	for (size_t i = 0; i < resolution; ++i) {
		double u = (double)i / (double)(resolution - 1);

		for (size_t j = 0; j < resolution; ++j) {
			double v = (double)j / (double)(resolution - 1);

			bezierRect.generateBernsteinDerivatives(u, v, bernsteinDerivativesU, bernsteinDerivativesV);

			handles.push_back(mesh.add_vertex(Vector(static_cast<double*>(bezierRect.position(bernsteinDerivativesU, bernsteinDerivativesV)))));
			means.push_back(bezierRect.meanCurvature(bernsteinDerivativesU, bernsteinDerivativesV));
		}
	}

	//Creating triangles
	for (size_t i = 0; i < resolution - 1; ++i)
		for (size_t j = 0; j < resolution - 1; ++j) {
			tri.clear();
			tri.push_back(handles[i * resolution + j]);
			tri.push_back(handles[i * resolution + j + 1]);
			tri.push_back(handles[(i + 1) * resolution + j]);
			mesh.add_face(tri);
			tri.clear();
			tri.push_back(handles[(i + 1) * resolution + j]);
			tri.push_back(handles[i * resolution + j + 1]);
			tri.push_back(handles[(i + 1) * resolution + j + 1]);
			mesh.add_face(tri);
		}

	//Saving Mean Curvature to vertexes
	for (size_t i = 0; i < resolution * resolution; i++) {
		mesh.data(handles[i]).mean = means[i];
	}

}

//Generate Bezier Triangle Mesh and Mean curvature
void MyViewer::generateBezierTriMesh(int index, bool add) {
	size_t resolution = 50;
	if (!add) {
		mesh.clear();
	}
	std::vector<MyMesh::VertexHandle> tri;
	std::vector < std::vector<MyMesh::VertexHandle> > handles;
	std::vector < std::vector< double > > means;

	std::vector< std::vector<std::vector<double>> > bernsteinDerivatives = std::vector < std::vector<std::vector<double>>>();

	//Calculate vertexes and their curvature
	for (size_t i = 0; i <= resolution; ++i) {
		handles.push_back(std::vector<MyMesh::VertexHandle>());
		means.push_back(std::vector<double>());
		for (size_t j = 0; j <= i; ++j) {
			double u = (double)(i - j) / (double)(resolution);
			double v = (double)j / (double)(resolution);

			bezierTri[index].generateBernsteinDerivatives(u, v, bernsteinDerivatives);

			handles[i].push_back(mesh.add_vertex(Vector(static_cast<double*>(bezierTri[index].position(bernsteinDerivatives)))));
			means[i].push_back(bezierTri[index].meanCurvature(bernsteinDerivatives));
		}
	}

	//Creating triangles
	for (size_t i = 0; i < resolution; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			tri.clear();
			tri.push_back(handles[i][j]);
			tri.push_back(handles[i + 1][j + 1]);
			tri.push_back(handles[i + 1][j]);
			mesh.add_face(tri);
			if (i != resolution - 1) {
				tri.clear();
				tri.push_back(handles[i + 1][j]);
				tri.push_back(handles[i + 1][j + 1]);
				tri.push_back(handles[i + 2][j + 1]);
				mesh.add_face(tri);
			}
		}
	}

	//Saving Mean Curvature to vertexes
	for (size_t i = 0; i <= resolution; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			mesh.data(handles[i][j]).mean = means[i][j] * (1 - 2 * index);
			//mesh.data(handles[i][j]).mean = means[i][(index*i)+j*(1-2*index)] * (1 - 2 * index);
		}
	}
}
void MyViewer::generateBezierTriPairMesh() {
	generateBezierTriMesh(Master, false);
	generateBezierTriMesh(Slave, true);
}

//Elevate degree of Bezier Rectangle
void MyViewer::elevateBezierRectDegree() {
	bezierRect.elevateDegree();
}
//Elevate degree of Bezier Triangle
void MyViewer::elevateBezierTriDegree(int index) {
	bezierTri[index].elevateDegree();
}
//Elevate degree of Bezier Triangle
void MyViewer::elevateBezierTriPairDegree() {
	bezierTri[Master].elevateDegree();
	bezierTri[Slave].elevateDegree();
}

//Drawing control net of Bezier Rectangle
void MyViewer::drawBezierRectControlNet() const {
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);
	glColor3d(0.3, 0.3, 1.0);

	size_t degree[2];
	degree[0] = bezierRect.getK();
	degree[1] = bezierRect.getL();
	for (size_t chooser = 0; chooser < 2; ++chooser)
		for (size_t i = 0; i <= degree[chooser]; ++i) {
			glBegin(GL_LINE_STRIP);
			for (size_t j = 0; j <= degree[1 - chooser]; ++j) {
				size_t const k = chooser ? j : i;
				size_t const l = chooser ? i : j;
				const auto& p = bezierRect.getControlPoint(k, l);
				glVertex3dv(p);
			}
			glEnd();
		}
	glLineWidth(1.0);
	glPointSize(8.0);
	glColor3d(1.0, 0.0, 1.0);
	glBegin(GL_POINTS);
	for (size_t k = 0; k <= bezierRect.getK(); ++k)
		for (size_t l = 0; l <= bezierRect.getL(); ++l)
			glVertex3dv(bezierRect.getControlPoint(k, l));
	glEnd();
	glPointSize(1.0);
	glEnable(GL_LIGHTING);
}

//Drawing control net of Bezier Triangle
void MyViewer::drawBezierTriControlNet(int index) const {
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);
	glColor3d(0.3, 0.3 + index / 2.0, 1.0 - index);

	//(1,0,0)  direction
	for (size_t i = 0; i <= bezierTri[index].getDegree(); ++i) {
		glBegin(GL_LINE_STRIP);
		for (size_t j = 0; j <= i; ++j) {
			int row = i;
			int column = j;
			int k = row - column;
			int l = column;
			const auto& p = bezierTri[index].getControlPoint(k, l);
			glVertex3dv(p);
		}
		glEnd();
	}
	//(0,0,1)  direction
	for (size_t i = 0; i <= bezierTri[index].getDegree(); ++i) {
		glBegin(GL_LINE_STRIP);
		for (size_t j = 0; j <= i; ++j) {
			int row = bezierTri[index].getDegree() - j;
			int column = bezierTri[index].getDegree() - i;
			int k = row - column;
			int l = column;
			const auto& p = bezierTri[index].getControlPoint(k, l);
			glVertex3dv(p);
		}
		glEnd();
	}
	//(0,1,0)  direction
	for (size_t i = 0; i <= bezierTri[index].getDegree(); ++i) {
		glBegin(GL_LINE_STRIP);
		for (size_t j = 0; j <= i; ++j) {
			int row = bezierTri[index].getDegree() - i + j;
			int column = j;
			int k = row - column;
			int l = column;
			const auto& p = bezierTri[index].getControlPoint(k, l);
			glVertex3dv(p);
		}
		glEnd();
	}

	glLineWidth(1.0);
	glPointSize(8.0);
	glColor3d(1.0, 0.0, 1.0);
	glBegin(GL_POINTS);
	for (size_t row = 0; row <= bezierTri[index].getDegree(); row++)
		for (size_t column = 0; column <= row; column++) {
			int k = row - column;
			int l = column;
			glVertex3dv(bezierTri[index].getControlPoint(k, l));
		}
	glEnd();
	glPointSize(1.0);
	glEnable(GL_LIGHTING);
}
void MyViewer::drawBezierTriPairControlNet() const {
	drawBezierTriControlNet(Master);
	drawBezierTriControlNet(Slave);
}

void MyViewer::updateMeanMinMax() {
	size_t n = mesh.n_vertices();
	if (n == 0)
		return;

	std::vector<double> mean;
	mean.reserve(n);
	for (auto v : mesh.vertices())
		mean.push_back(mesh.data(v).mean);

	std::sort(mean.begin(), mean.end());
	size_t k = (double)n * cutoff_ratio;
	mean_min = std::min(mean[k ? k - 1 : 0], 0.0);
	mean_max = std::max(mean[n - k], 0.0);
}

void MyViewer::localSystem(const MyViewer::Vector& normal,
	MyViewer::Vector& u, MyViewer::Vector& v) {
	// Generates an orthogonal (u,v) coordinate system in the plane defined by `normal`.
	int maxi = 0, nexti = 1;
	double max = std::abs(normal[0]), next = std::abs(normal[1]);
	if (max < next) {
		std::swap(max, next);
		std::swap(maxi, nexti);
	}
	if (std::abs(normal[2]) > max) {
		nexti = maxi;
		maxi = 2;
	}
	else if (std::abs(normal[2]) > next)
		nexti = 2;

	u.vectorize(0.0);
	u[nexti] = -normal[maxi];
	u[maxi] = normal[nexti];
	u /= u.norm();
	v = normal % u;
}

double MyViewer::voronoiWeight(MyViewer::MyMesh::HalfedgeHandle in_he) {
	// Returns the area of the triangle bounded by in_he that is closest
	// to the vertex pointed to by in_he.
	if (mesh.is_boundary(in_he))
		return 0;
	auto next = mesh.next_halfedge_handle(in_he);
	auto prev = mesh.prev_halfedge_handle(in_he);
	double c2 = mesh.calc_edge_vector(in_he).sqrnorm();
	double b2 = mesh.calc_edge_vector(next).sqrnorm();
	double a2 = mesh.calc_edge_vector(prev).sqrnorm();
	double alpha = mesh.calc_sector_angle(in_he);

	if (a2 + b2 < c2)                // obtuse gamma
		return 0.125 * b2 * std::tan(alpha);
	if (a2 + c2 < b2)                // obtuse beta
		return 0.125 * c2 * std::tan(alpha);
	if (b2 + c2 < a2) {              // obtuse alpha
		double b = std::sqrt(b2), c = std::sqrt(c2);
		double total_area = 0.5 * b * c * std::sin(alpha);
		double beta = mesh.calc_sector_angle(prev);
		double gamma = mesh.calc_sector_angle(next);
		return total_area - 0.125 * (b2 * std::tan(gamma) + c2 * std::tan(beta));
	}

	double r2 = 0.25 * a2 / std::pow(std::sin(alpha), 2); // squared circumradius
	auto area = [r2](double x2) {
		return 0.125 * std::sqrt(x2) * std::sqrt(std::max(4.0 * r2 - x2, 0.0));
	};
	return area(b2) + area(c2);
}

#ifndef BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
	if (model_type != ModelType::BEZIER_SURFACE_RECTANGULAR
		&& model_type != ModelType::BEZIER_SURFACE_TRIANGULAR
		&& model_type != ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR) {
		std::map<MyMesh::FaceHandle, double> face_area;
		std::map<MyMesh::VertexHandle, double> vertex_area;

		for (auto f : mesh.faces())
			face_area[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

		// Compute triangle strip areas
		for (auto v : mesh.vertices()) {
			vertex_area[v] = 0;
			mesh.data(v).mean = 0;
			for (auto f : mesh.vf_range(v))
				vertex_area[v] += face_area[f];
			vertex_area[v] /= 3.0;
		}

		// Compute mean values using dihedral angles
		for (auto v : mesh.vertices()) {
			for (auto h : mesh.vih_range(v)) {
				auto vec = mesh.calc_edge_vector(h);
				double angle = mesh.calc_dihedral_angle(h); // signed; returns 0 at the boundary
				mesh.data(v).mean += angle * vec.norm();
			}
			mesh.data(v).mean *= 0.25 / vertex_area[v];
		}
	}
	if (update_min_max)
		updateMeanMinMax();
}
#else // BETTER_MEAN_CURVATURE
void MyViewer::updateMeanCurvature(bool update_min_max) {
	// As in the paper:
	//   S. Rusinkiewicz, Estimating curvatures and their derivatives on triangle meshes.
	//     3D Data Processing, Visualization and Transmission, IEEE, 2004.

	std::map<MyMesh::VertexHandle, Vector> efgp; // 2nd principal form
	std::map<MyMesh::VertexHandle, double> wp;   // accumulated weight

	// Initial setup
	for (auto v : mesh.vertices()) {
		efgp[v].vectorize(0.0);
		wp[v] = 0.0;
	}

	for (auto f : mesh.faces()) {
		// Setup local edges, vertices and normals
		auto h0 = mesh.halfedge_handle(f);
		auto h1 = mesh.next_halfedge_handle(h0);
		auto h2 = mesh.next_halfedge_handle(h1);
		auto e0 = mesh.calc_edge_vector(h0);
		auto e1 = mesh.calc_edge_vector(h1);
		auto e2 = mesh.calc_edge_vector(h2);
		auto n0 = mesh.normal(mesh.to_vertex_handle(h1));
		auto n1 = mesh.normal(mesh.to_vertex_handle(h2));
		auto n2 = mesh.normal(mesh.to_vertex_handle(h0));

		Vector n = mesh.normal(f), u, v;
		localSystem(n, u, v);

		// Solve a LSQ equation for (e,f,g) of the face
		Eigen::MatrixXd A(6, 3);
		A << (e0 | u), (e0 | v), 0.0,
			0.0, (e0 | u), (e0 | v),
			(e1 | u), (e1 | v), 0.0,
			0.0, (e1 | u), (e1 | v),
			(e2 | u), (e2 | v), 0.0,
			0.0, (e2 | u), (e2 | v);
		Eigen::VectorXd b(6);
		b << ((n2 - n1) | u),
			((n2 - n1) | v),
			((n0 - n2) | u),
			((n0 - n2) | v),
			((n1 - n0) | u),
			((n1 - n0) | v);
		Eigen::Vector3d x = A.fullPivLu().solve(b);

		Eigen::Matrix2d F;          // Fundamental matrix for the face
		F << x(0), x(1),
			x(1), x(2);

		for (auto h : mesh.fh_range(f)) {
			auto p = mesh.to_vertex_handle(h);

			// Rotate the (up,vp) local coordinate system to be coplanar with that of the face
			Vector np = mesh.normal(p), up, vp;
			localSystem(np, up, vp);
			auto axis = (np % n).normalize();
			double angle = std::acos(std::min(std::max(n | np, -1.0), 1.0));
			auto rotation = Eigen::AngleAxisd(angle, Eigen::Vector3d(axis.data()));
			Eigen::Vector3d up1(up.data()), vp1(vp.data());
			up1 = rotation * up1;    vp1 = rotation * vp1;
			up = Vector(up1.data()); vp = Vector(vp1.data());

			// Compute the vertex-local (e,f,g)
			double e, f, g;
			Eigen::Vector2d upf, vpf;
			upf << (up | u), (up | v);
			vpf << (vp | u), (vp | v);
			e = upf.transpose() * F * upf;
			f = upf.transpose() * F * vpf;
			g = vpf.transpose() * F * vpf;

			// Accumulate the results with Voronoi weights
			double w = voronoiWeight(h);
			efgp[p] += Vector(e, f, g) * w;
			wp[p] += w;
		}
	}

	// Compute the principal curvatures
	for (auto v : mesh.vertices()) {
		auto& efg = efgp[v];
		efg /= wp[v];
		Eigen::Matrix2d F;
		F << efg[0], efg[1],
			efg[1], efg[2];
		auto k = F.eigenvalues();   // always real, because F is a symmetric real matrix
		mesh.data(v).mean = (k(0).real() + k(1).real()) / 2.0;
	}

	if (update_min_max)
		updateMeanMinMax();
}
#endif

static Vec HSV2RGB(Vec hsv) {
	// As in Wikipedia
	double c = hsv[2] * hsv[1];
	double h = hsv[0] / 60;
	double x = c * (1 - std::abs(std::fmod(h, 2) - 1));
	double m = hsv[2] - c;
	Vec rgb(m, m, m);
	if (h <= 1)
		return rgb + Vec(c, x, 0);
	if (h <= 2)
		return rgb + Vec(x, c, 0);
	if (h <= 3)
		return rgb + Vec(0, c, x);
	if (h <= 4)
		return rgb + Vec(0, x, c);
	if (h <= 5)
		return rgb + Vec(x, 0, c);
	if (h <= 6)
		return rgb + Vec(c, 0, x);
	return rgb;
}

Vec MyViewer::meanMapColor(double d) const {
	double red = 0, green = 120, blue = 240; // Hue
	if (d < 0) {
		double alpha = mean_min ? std::min(d / mean_min, 1.0) : 1.0;
		return HSV2RGB({ green * (1 - alpha) + blue * alpha, 1, 1 });
	}
	double alpha = mean_max ? std::min(d / mean_max, 1.0) : 1.0;
	return HSV2RGB({ green * (1 - alpha) + red * alpha, 1, 1 });
}

void MyViewer::fairMesh() {
	if (model_type != ModelType::MESH)
		return;

	emit startComputation(tr("Fairing mesh..."));
	OpenMesh::Smoother::JacobiLaplaceSmootherT<MyMesh> smoother(mesh);
	smoother.initialize(OpenMesh::Smoother::SmootherT<MyMesh>::Normal, // or: Tangential_and_Normal
		OpenMesh::Smoother::SmootherT<MyMesh>::C1);
	for (size_t i = 1; i <= 10; ++i) {
		smoother.smooth(10);
		emit midComputation(i * 10);
	}
	updateMesh(false);
	emit endComputation();
}

void MyViewer::updateVertexNormals() {
	// Weights according to:
	//   N. Max, Weights for computing vertex normals from facet normals.
	//     Journal of Graphics Tools, Vol. 4(2), 1999.
	for (auto v : mesh.vertices()) {
		Vector n(0.0, 0.0, 0.0);
		for (auto h : mesh.vih_range(v)) {
			if (mesh.is_boundary(h))
				continue;
			auto in_vec = mesh.calc_edge_vector(h);
			auto out_vec = mesh.calc_edge_vector(mesh.next_halfedge_handle(h));
			double w = in_vec.sqrnorm() * out_vec.sqrnorm();
			n += (in_vec % out_vec) / (w == 0.0 ? 1.0 : w);
		}
		double len = n.length();
		if (len != 0.0)
			n /= len;
		mesh.set_normal(v, n);
	}
}

void MyViewer::updateMesh(bool update_mean_range) {
	if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR)
		generateBezierRectMesh();
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR)
		generateBezierTriMesh();
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR)
		generateBezierTriPairMesh();
	if (visualization == Visualization::HOMEWORK)
		recalculateQuarterSizes();
	mesh.request_face_normals(); mesh.request_vertex_normals();
	mesh.update_face_normals(); //mesh.update_vertex_normals();
	updateVertexNormals();
	updateMeanCurvature(update_mean_range);
}

void MyViewer::setupCamera() {
	// Set camera on the model
	Vector box_min, box_max;
	box_min = box_max = mesh.point(*mesh.vertices_begin());
	for (auto v : mesh.vertices()) {
		box_min.minimize(mesh.point(v));
		box_max.maximize(mesh.point(v));
	}
	camera()->setSceneBoundingBox(Vec(box_min.data()), Vec(box_max.data()));
	camera()->showEntireScene();

	slicing_scaling = 20 / (box_max - box_min).max();

	setSelectedName(-1);
	axes.shown = false;

	update();
}

bool MyViewer::openMesh(const std::string& filename, bool update_view) {
	if (!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
		return false;
	model_type = ModelType::MESH;
	last_filename = filename;
	updateMesh(update_view);
	if (update_view)
		setupCamera();
	return true;
}

void MyViewer::init() {
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	QImage img(":/isophotes.png");
	glGenTextures(1, &isophote_texture);
	glBindTexture(GL_TEXTURE_2D, isophote_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img.width(), img.height(), 0, GL_BGRA,
		GL_UNSIGNED_BYTE, img.convertToFormat(QImage::Format_ARGB32).bits());

	QImage img2(":/environment.png");
	glGenTextures(1, &environment_texture);
	glBindTexture(GL_TEXTURE_2D, environment_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img2.width(), img2.height(), 0, GL_BGRA,
		GL_UNSIGNED_BYTE, img2.convertToFormat(QImage::Format_ARGB32).bits());

	glGenTextures(1, &slicing_texture);
	glBindTexture(GL_TEXTURE_1D, slicing_texture);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	static const unsigned char slicing_img[] = { 0b11111111, 0b00011100 };
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 2, 0, GL_RGB, GL_UNSIGNED_BYTE_3_3_2, &slicing_img);
}

void MyViewer::recalculateQuarterSizes() {
	auto areas = std::vector<double>();

	for (auto f : mesh.faces()) {
		mesh.data(f).area =
			mesh.calc_sector_area(
				mesh.fh_range(f).begin()
			);
		areas.push_back(mesh.data(f).area);
	}

	double tmp;
	std::sort(areas.begin(), areas.end(), std::greater<double>());

	biggerQuarter = areas[areas.size() / 4];
	emit median(sqrt(areas[areas.size() / 2]));
	smallerQuarter = areas[areas.size() * 3 / 4];
}



void MyViewer::draw() {
	if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR && show_control_points)
		drawBezierRectControlNet();
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR && show_control_points)
		drawBezierTriControlNet();
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR && show_control_points)
		drawBezierTriPairControlNet();

	glPolygonMode(GL_FRONT_AND_BACK, !show_solid && show_wireframe ? GL_LINE : GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1, 1);

	if (show_solid || show_wireframe) {
		if (visualization == Visualization::PLAIN)
			glColor3d(1.0, 1.0, 1.0);
		else if (visualization == Visualization::ISOPHOTES) {
			glBindTexture(GL_TEXTURE_2D, current_isophote_texture);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glEnable(GL_TEXTURE_2D);
			glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
			glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_SPHERE_MAP);
			glEnable(GL_TEXTURE_GEN_S);
			glEnable(GL_TEXTURE_GEN_T);
		}
		else if (visualization == Visualization::SLICING) {
			glBindTexture(GL_TEXTURE_1D, slicing_texture);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glEnable(GL_TEXTURE_1D);
		}

		for (auto f : mesh.faces()) {

			if (visualization == Visualization::HOMEWORK) {

				double area = mesh.data(f).area;
				if (area <= biggerQuarter && area > smallerQuarter) {
					glColor3d(1.0, 0.0, 0.0);
				}
				else {
					glColor3d(1.0, 1.0, 1.0);
				}
			}
			else {
				glColor3d(1.0, 1.0, 1.0);
			}

			glBegin(GL_POLYGON);
			for (auto v : mesh.fv_range(f)) {
				if (visualization == Visualization::MEAN)
					glColor3dv(meanMapColor(mesh.data(v).mean));
				else if (visualization == Visualization::SLICING)
					glTexCoord1d(mesh.point(v) | slicing_dir * slicing_scaling);
				glNormal3dv(mesh.normal(v).data());
				glVertex3dv(mesh.point(v).data());
			}
			glEnd();
		}
		if (visualization == Visualization::ISOPHOTES) {
			glDisable(GL_TEXTURE_GEN_S);
			glDisable(GL_TEXTURE_GEN_T);
			glDisable(GL_TEXTURE_2D);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		}
		else if (visualization == Visualization::SLICING) {
			glDisable(GL_TEXTURE_1D);
		}
	}

	if (show_solid && show_wireframe) {
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3d(0.0, 0.0, 0.0);
		glDisable(GL_LIGHTING);
		for (auto f : mesh.faces()) {
			glBegin(GL_POLYGON);
			for (auto v : mesh.fv_range(f))
				glVertex3dv(mesh.point(v).data());
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}

	if (axes.shown)
		drawAxes();
}

void MyViewer::drawAxes() const {
	const Vec& p = axes.position;
	glColor3d(1.0, 0.0, 0.0);
	drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
	glColor3d(0.0, 1.0, 0.0);
	drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
	glColor3d(0.0, 0.0, 1.0);
	drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
	glEnd();
}

void MyViewer::drawWithNames() {
	if (axes.shown)
		return drawAxesWithNames();

	switch (model_type) {
	case ModelType::NONE: break;
	case ModelType::MESH:
		if (!show_wireframe)
			return;
		for (auto v : mesh.vertices()) {
			glPushName(v.idx());
			glRasterPos3dv(mesh.point(v).data());
			glPopName();
		}
		break;
	case ModelType::BEZIER_SURFACE_RECTANGULAR:
		if (!show_control_points)
			return;
		for (size_t k = 0; k <= bezierRect.getK(); k++) {
			for (size_t l = 0; l <= bezierRect.getL(); l++) {
				Vec const& p = bezierRect.getControlPoint(k, l);
				glPushName(k * (bezierRect.getL() + 1) + l);
				glRasterPos3fv(p);
				glPopName();
			}
		}
		break;
	case ModelType::BEZIER_SURFACE_TRIANGULAR:
		if (!show_control_points)
			return;
		for (size_t row = 0; row <= bezierTri[Master].getDegree(); row++)
			for (size_t column = 0; column <= row; column++) {
				int k = row - column;
				int l = column;
				Vec const& p = bezierTri[Master].getControlPoint(k, l);
				glPushName(row * (row + 1) / 2 + column);
				glRasterPos3fv(p);
				glPopName();
			}
	case ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR:
		if (!show_control_points)
			return;
		for (int index = 0; index < 2; index++)
			for (size_t row = 0; row <= bezierTri[index].getDegree(); row++)
				for (size_t column = 0; column <= row; column++) {
					int k = row - column;
					int l = column;
					Vec const& p = bezierTri[index].getControlPoint(k, l);
					glPushName(
						row * (row + 1) / 2
						+ column
						+ (index
							* (bezierTri[Master].getDegree() + 1)
							* (bezierTri[Master].getDegree() + 2)
							/ 2)
					);
					glRasterPos3fv(p);
					glPopName();
				}
		break;
	}
}

void MyViewer::drawAxesWithNames() const {
	const Vec& p = axes.position;
	glPushName(0);
	drawArrow(p, p + Vec(axes.size, 0.0, 0.0), axes.size / 50.0);
	glPopName();
	glPushName(1);
	drawArrow(p, p + Vec(0.0, axes.size, 0.0), axes.size / 50.0);
	glPopName();
	glPushName(2);
	drawArrow(p, p + Vec(0.0, 0.0, axes.size), axes.size / 50.0);
	glPopName();
}

void MyViewer::postSelection(const QPoint& p) {
	int sel = selectedName();
	if (sel == -1) {
		axes.shown = false;
		return;
	}

	if (axes.shown) {
		axes.selected_axis = sel;
		bool found;
		axes.grabbed_pos = camera()->pointUnderPixel(p, found);
		axes.original_pos = axes.position;
		if (!found)
			axes.shown = false;
		return;
	}

	selected_vertex = sel;
	if (model_type == ModelType::MESH)
		axes.position = Vec(mesh.point(MyMesh::VertexHandle(sel)).data());
	if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR)
		axes.position = bezierRect.getControlPoint(sel);
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR)
		axes.position = bezierTri[Master].getControlPoint(sel);
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR) {
		int masterSize = (bezierTri[Master].getDegree() + 1) * (bezierTri[Master].getDegree() + 2) / 2;
		int index = sel / masterSize;
		int kmi = sel;
		if (index == Slave)
			kmi -= masterSize;
		axes.position = bezierTri[index].getControlPoint(kmi);
	}
	double depth = camera()->projectedCoordinatesOf(axes.position)[2];
	Vec q1 = camera()->unprojectedCoordinatesOf(Vec(0.0, 0.0, depth));
	Vec q2 = camera()->unprojectedCoordinatesOf(Vec(width(), height(), depth));
	axes.size = (q1 - q2).norm() / 10.0;
	axes.shown = true;
	axes.selected_axis = -1;
}

void MyViewer::keyPressEvent(QKeyEvent* e) {
	OpenMesh::Subdivider::Uniform::Sqrt3T<MyMesh> sqrtThree;
	OpenMesh::Subdivider::Uniform::LoopT<MyMesh> loop;
	if (e->modifiers() == Qt::NoModifier)
		switch (e->key()) {
		case Qt::Key_R:
			if (model_type == ModelType::MESH)
				openMesh(last_filename, false);
			else if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR)
				openBezierRect(last_filename, false);
			else if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR)
				openBezierTri(last_filename, false);
			else if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR)
				openBezierTriPair(last_filename, false);
			update();
			break;
		case Qt::Key_O:
			if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
				camera()->setType(qglviewer::Camera::ORTHOGRAPHIC);
			else
				camera()->setType(qglviewer::Camera::PERSPECTIVE);
			update();
			break;
		case Qt::Key_P:
			visualization = Visualization::PLAIN;
			update();
			break;
		case Qt::Key_M:
			visualization = Visualization::MEAN;
			update();
			break;
		case Qt::Key_L:
			visualization = Visualization::SLICING;
			update();
			break;
		case Qt::Key_I:
			visualization = Visualization::ISOPHOTES;
			current_isophote_texture = isophote_texture;
			update();
			break;
		case Qt::Key_C:
			show_control_points = !show_control_points;
			update();
			break;
		case Qt::Key_S:
			show_solid = !show_solid;
			update();
			break;
		case Qt::Key_W:
			show_wireframe = !show_wireframe;
			update();
			break;
		case Qt::Key_F:
			fairMesh();
			update();
			break;
		case Qt::Key_G:
			continuity = !continuity;
			break;
		case Qt::Key_E:
			if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR) {
				elevateBezierRectDegree();
				updateMesh(true);
				update();
			}
			else if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR) {
				elevateBezierTriDegree();
				updateMesh(true);
				update();
			}
			else if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR) {
				elevateBezierTriPairDegree();
				updateMesh(true);
				update();
			}
			break;
			/*
			case Qt::Key_X:
			  visualization = Visualization::HOMEWORK;
			  updateMesh();
			  update();
			  break;
			  */
		case Qt::Key_3:
			// Execute subdivision step

			sqrtThree.attach(mesh);
			sqrtThree(1);
			sqrtThree.detach();
			updateMesh(true);
			update();
			break;
		default:
			QGLViewer::keyPressEvent(e);
		}
	else if (e->modifiers() == Qt::KeypadModifier)
		switch (e->key()) {
		case Qt::Key_Plus:
			slicing_scaling *= 2;
			update();
			break;
		case Qt::Key_Minus:
			slicing_scaling /= 2;
			update();
			break;
		case Qt::Key_Asterisk:
			slicing_dir = Vector(static_cast<double*>(camera()->viewDirection()));
			update();
			break;
		}
	else if (e->modifiers() == Qt::ShiftModifier) {
		// Initialize subdivider
		switch (e->key()) {
		case Qt::Key_L:
			// Execute subdivision step
			loop.attach(mesh);
			loop(1);
			loop.detach();
			updateMesh(true);
			update();
			break;
		default:
			QGLViewer::keyPressEvent(e);
		}
	}
	else
		QGLViewer::keyPressEvent(e);
}

Vec MyViewer::intersectLines(const Vec& ap, const Vec& ad, const Vec& bp, const Vec& bd) {
	// always returns a point on the (ap, ad) line
	double a = ad * ad, b = ad * bd, c = bd * bd;
	double d = ad * (ap - bp), e = bd * (ap - bp);
	if (a * c - b * b < 1.0e-7)
		return ap;
	double s = (b * e - c * d) / (a * c - b * b);
	return ap + s * ad;
}

void MyViewer::mouseMoveEvent(QMouseEvent* e) {
	if (!axes.shown ||
		(axes.selected_axis < 0 && !(e->modifiers() & Qt::ControlModifier)) ||
		!(e->modifiers() & (Qt::ShiftModifier | Qt::ControlModifier)) ||
		!(e->buttons() & Qt::LeftButton))
		return QGLViewer::mouseMoveEvent(e);

	if (e->modifiers() & Qt::ControlModifier) {
		// move in screen plane
		double depth = camera()->projectedCoordinatesOf(axes.position)[2];
		axes.position = camera()->unprojectedCoordinatesOf(Vec(e->pos().x(), e->pos().y(), depth));
	}
	else {
		Vec from, dir, axis(axes.selected_axis == 0, axes.selected_axis == 1, axes.selected_axis == 2);
		camera()->convertClickToLine(e->pos(), from, dir);
		auto p = intersectLines(axes.grabbed_pos, axis, from, dir);
		float d = (p - axes.grabbed_pos) * axis;
		axes.position[axes.selected_axis] = axes.original_pos[axes.selected_axis] + d;
	}

	if (model_type == ModelType::MESH)
		mesh.set_point(MyMesh::VertexHandle(selected_vertex),
			Vector(static_cast<double*>(axes.position)));
	if (model_type == ModelType::BEZIER_SURFACE_RECTANGULAR)
		bezierRect.setControlPoint(selected_vertex, axes.position);
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR)
		bezierTri[Master].setControlPoint(selected_vertex, axes.position);
	if (model_type == ModelType::BEZIER_SURFACE_TRIANGULAR_PAIR) {
		int masterSize = (bezierTri[Master].getDegree() + 1) * (bezierTri[Master].getDegree() + 2) / 2;
		int index = selected_vertex / masterSize;
		int kmi = selected_vertex;
		if (index == Slave)
			kmi -= masterSize;
		Vec alfaBeta0;
		Vec alfaBeta1;
		if (continuity) {
			alfaBeta0 = bezierTri[Slave].updateAlfaBeta(
				bezierTri[Slave].getControlPoint(bezierTri[Slave].getDegree(), 0),
				bezierTri[Slave].getControlPoint(bezierTri[Slave].getDegree() - 1, 1),
				bezierTri[Master].getControlPoint(0, bezierTri[Master].getDegree() - 1),
				bezierTri[Slave].getControlPoint(bezierTri[Slave].getDegree() - 1, 0)
			);
			alfaBeta1 = bezierTri[Slave].updateAlfaBeta(
				bezierTri[Slave].getControlPoint(1, bezierTri[Slave].getDegree() - 1),
				bezierTri[Slave].getControlPoint(0, bezierTri[Slave].getDegree()),
				bezierTri[Master].getControlPoint(bezierTri[Slave].getDegree() - 1, 0),
				bezierTri[Slave].getControlPoint(0, bezierTri[Master].getDegree() - 1)
			);
		}
		bezierTri[index].setControlPoint(kmi, axes.position);
		if (continuity)
			bezierTri[Slave].updateSlave(bezierTri[Master], alfaBeta0[0], alfaBeta1[0], alfaBeta0[1], alfaBeta1[1]);
	}
	updateMesh();
	update();
}


QString MyViewer::helpString() const {
	QString text("<h2>Sample Framework</h2>"
		"<p>This is a minimal framework for 3D mesh manipulation, which can be "
		"extended and used as a base for various projects, for example "
		"prototypes for fairing algorithms, or even displaying/modifying "
		"parametric surfaces, etc.</p>"
		"<p>The following hotkeys are available:</p>"
		"<ul>"
		"<li>&nbsp;R: Reload model</li>"
		"<li>&nbsp;O: Toggle orthographic projection</li>"
		"<li>&nbsp;P: Set plain map (no coloring)</li>"
		"<li>&nbsp;M: Set mean curvature map</li>"
		"<li>&nbsp;L: Set slicing map<ul>"
		"<li>&nbsp;+: Increase slicing density</li>"
		"<li>&nbsp;-: Decrease slicing density</li>"
		"<li>&nbsp;*: Set slicing direction to view</li></ul></li>"
		"<li>&nbsp;I: Set isophote line map</li>"
		"<li>&nbsp;E: Set environment texture</li>"
		"<li>&nbsp;C: Toggle control polygon visualization</li>"
		"<li>&nbsp;S: Toggle solid (filled polygon) visualization</li>"
		"<li>&nbsp;W: Toggle wireframe visualization</li>"
		"<li>&nbsp;F: Fair mesh</li>"
		"</ul>"
		"<p>There is also a simple selection and movement interface, enabled "
		"only when the wireframe/controlnet is displayed: a mesh vertex can be selected "
		"by shift-clicking, and it can be moved by shift-dragging one of the "
		"displayed axes. Pressing ctrl enables movement in the screen plane.</p>"
		"<p>Note that libQGLViewer is furnished with a lot of useful features, "
		"such as storing/loading view positions, or saving screenshots. "
		"OpenMesh also has a nice collection of tools for mesh manipulation: "
		"decimation, subdivision, smoothing, etc. These can provide "
		"good comparisons to the methods you implement.</p>"
		"<p>This software can be used as a sample GUI base for handling "
		"parametric or procedural surfaces, as well. The power of "
		"Qt and libQGLViewer makes it easy to set up a prototype application. "
		"Feel free to modify and explore!</p>"
		"<p align=\"right\">Peter Salvi</p>");
	return text;
}
