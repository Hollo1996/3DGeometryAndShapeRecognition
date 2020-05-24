// -*- mode: c++ -*-
#pragma once

#include <string>

#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/LoopT.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/Sqrt3T.hh>
#include <BezierTri.h>
#include <BezierRect.h>

using qglviewer::Vec;

using qglviewer::Vec;

class MyViewer : public QGLViewer {
  Q_OBJECT

  double biggerQuarter;
  double smallerQuarter;

public:
  explicit MyViewer(QWidget *parent);
  virtual ~MyViewer();

  inline double getCutoffRatio() const;
  inline void setCutoffRatio(double ratio);

  inline double getMeanMin() const;
  inline void setMeanMin(double min);
  inline double getMeanMax() const;
  inline void setMeanMax(double max);

  inline const double *getSlicingDir() const;
  inline void setSlicingDir(double x, double y, double z);
  inline double getSlicingScaling() const;
  inline void setSlicingScaling(double scaling);

  bool openMesh(const std::string &filename, bool update_view = true);

  //BEZIER
  bool openBezierRect(const std::string& filename, bool update_view = true);
  bool saveBezierRect(const std::string& filename);
  bool openBezierTri(const std::string& filename, bool update_view = true);
  bool saveBezierTri(const std::string& filename);
  bool openBezierTriPair(const std::string& filename, bool update_view = true);
  bool saveBezierTriPair(const std::string& filename);

signals:
  void startComputation(QString message);
  void midComputation(int percent);
  void endComputation();
  void median(double median);

protected:
  virtual void init() override;
  virtual void draw() override;
  virtual void drawWithNames() override;
  virtual void postSelection(const QPoint &p) override;
  virtual void keyPressEvent(QKeyEvent *e) override;
  virtual void mouseMoveEvent(QMouseEvent *e) override;
  virtual QString helpString() const override;

private:
  struct MyTraits : public OpenMesh::DefaultTraits {
    using Point  = OpenMesh::Vec3d; // the default would be Vec3f
    using Normal = OpenMesh::Vec3d;
    VertexTraits {
      double mean;              // approximated mean curvature
    };
    FaceTraits{
        double area;
    };
  };
  using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
  using Vector = OpenMesh::VectorT<double,3>;
  OpenMesh::Subdivider::Uniform::Sqrt3T<MyMesh> sqrtThree;
  OpenMesh::Subdivider::Uniform::LoopT<MyMesh> loop;

  // Mesh
  void recalculateQuarterSizes();
  void updateMesh(bool update_mean_range = true);
  void updateVertexNormals();
  void localSystem(const Vector &normal, Vector &u, Vector &v);
  double voronoiWeight(MyMesh::HalfedgeHandle in_he);
  void updateMeanMinMax();
  void updateMeanCurvature(bool update_min_max = true);

  //BEZIER
  bool continuity = false;
  BezierRect bezierRect;
  void generateBezierRectMesh();
  void elevateBezierRectDegree();
  void drawBezierRectControlNet() const;
  BezierTri bezierTri[2];
  const int Master = 0;
  const int Slave = 1;
  void generateBezierTriMesh(int index = 0, bool add = false);
  void elevateBezierTriDegree(int index = 0);
  void drawBezierTriControlNet(int index = 0) const;
  void generateBezierTriPairMesh();
  void elevateBezierTriPairDegree();
  void drawBezierTriPairControlNet() const;

  // Visualization
  void setupCamera();
  Vec meanMapColor(double d) const;
  void drawAxes() const;
  void drawAxesWithNames() const;
  static Vec intersectLines(const Vec &ap, const Vec &ad, const Vec &bp, const Vec &bd);

  // Other
  void fairMesh();

  //////////////////////
  // Member variables //
  //////////////////////

  enum class ModelType { 
      NONE, 
      MESH, 
      BEZIER_SURFACE_RECTANGULAR, 
      BEZIER_SURFACE_TRIANGULAR, 
      BEZIER_SURFACE_TRIANGULAR_PAIR
  } model_type;

  // Mesh
  MyMesh mesh;
  size_t resolution;

  // Visualization
  double mean_min, mean_max, cutoff_ratio;
  bool show_control_points, show_solid, show_wireframe;
  enum class Visualization { PLAIN, MEAN, SLICING, ISOPHOTES, HOMEWORK } visualization;
  GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
  Vector slicing_dir;
  double slicing_scaling;
  int selected_vertex;
  struct ModificationAxes {
    bool shown;
    float size;
    int selected_axis;
    Vec position, grabbed_pos, original_pos;
  } axes;
  std::string last_filename;

};

#include "MyViewer.hpp"
