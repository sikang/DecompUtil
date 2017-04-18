#ifndef ELLIPSE_DECOMP_H
#define ELLIPSE_DECOMP_H

#define USE_CENTROID 1

#include <decomp_util/data_type.h>
#include <time.h>
#include <memory>
#include "line_segment.h"
#include "max_ellipsoid.hpp"

class EllipseDecomp {
public:
  EllipseDecomp() {};
  EllipseDecomp(const Vec3f &origin, const Vec3f &dim, bool verbose = false);
	~EllipseDecomp(){};
  void set_obstacles(const vec_Vec3f &obs) { obs_ = obs; }

  vec_Vec3f get_dilate_path() const { return dilate_path_; }
  vec_Vec3f get_center_path() const { return center_path_; }
  Polyhedra get_polyhedra(int id = 0) const;
  vec_Ellipsoid get_ellipsoids() const { return ellipsoids_; }
  vec_Vec3f get_pts() const;
  vec_LinearConstraint3f get_constraints();
  decimal_t get_corridor_volume();
  decimal_t get_ellipsoid_volume();

  void clean();
  void info();

  void change_end(decimal_t z_up = 0.1, decimal_t z_down = 0.1);
  bool decomp(const vec_Vec3f &poses,
      decimal_t d = 0.0,
      decimal_t ds = 1.0,
      decimal_t h = 3.0);

protected:
  void clear();
  void add_bounding(Polyhedron &Vs);

  vec_Vec3f cal_centers(const Polyhedra &Vs);

  vec_Vec3f dilate_path_;
  vec_Vec3f center_path_;
  vec_Vec3f obs_;

  vec_Ellipsoid ellipsoids_;
  Polyhedra polyhedrons_;
  Polyhedra intersect_polyhedrons_;
  std::vector<std::shared_ptr<LineSegment>> lines_;

  Vec3f min_; // bounding box for visualization
  Vec3f max_;

  bool verbose_;
};
#endif
