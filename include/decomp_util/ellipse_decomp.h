/**
 * @file ellipse_decomp.h
 * @brief EllipseDecomp Class
 */
#ifndef ELLIPSE_DECOMP_H
#define ELLIPSE_DECOMP_H

#define USE_CENTROID 1

#include <decomp_util/data_type.h>
#include <time.h>
#include <memory>
#include "line_segment.h"
#include "max_ellipsoid.hpp"

/**
 * @brief EllipseDecomp Class
 *
 * EllipseDecomp takes input as a given path and find the Safe Flight Corridor around it using Ellipsoids
 */
template <int Dim>
class EllipsoidDecomp {
public:
  ///Simple constructor
  EllipsoidDecomp(bool verbose = false);
  /**
   * @brief Basic constructor
   * @param origin The origin of the global bounding box
   * @param dim The dimension of the global bounding box
   */
  EllipsoidDecomp(const Vecf<Dim> &origin, const Vecf<Dim> &dim, bool verbose = false);
  ///Set obstacle points
  void set_obstacles(const vec_Vecf<Dim> &obs) { obs_ = obs; }
  ///Set dimension of bounding box
  void set_local_bounding_box(const Vecf<Dim>& bbox) { local_bbox_ = bbox; }

  ///Get the path that is used for dilation
  vec_Vecf<Dim> get_dilate_path() const { return dilate_path_; }
  ///Get the new path in the center of the Safe Flight Corridor
  vec_Vecf<Dim> get_center_path() const { return center_path_; }
  ///Get the Safe Flight Corridor
  Polyhedra<Dim> get_polyhedra() const { return polyheda_; }
  ///Get the intersected part of Safe Flight Corridor
  Polyhedra<Dim> get_intersect_polyhedra() const { return intersect_polyheda_; }
  ///Get the ellipsoids
  vec_E<Ellipsoid<Dim>> get_ellipsoids() const { return ellipsoids_; }
  ///Get the constraints of SFC as \f$Ax\leq b \f$
  vec_E<LinearConstraint<Dim>> get_constraints();
  ///Calculate the total volume of the SFC
  decimal_t get_corridor_volume();
  ///Calculate the total volume of the ellipsoids
  decimal_t get_ellipsoid_volume();

  ///Clean both dilation and obstacles, reset the whole class
  void clean();

  /**
   * @brief Decomposition thread
   * @param ppath The path to dilate
   * @param offset_x offset added to the long semi-axis, default is 0
   */
  bool decomp(const vec_Vecf<Dim> &path, double offset_x = 0);

  /**
   * @brief Shrink the safe flight corridor
   * @param path path that used to shirnk
   * @param shrink_distance shrinking distance, should be positive
   */
  void shrink(const vec_Vecf<Dim> &path, double shrink_distance);

protected:
  void clear();
  void add_bounding(Polyhedron<Dim> &Vs);

  vec_Vecf<Dim> cal_centers(const Polyhedra<Dim> &Vs);

  vec_Vecf<Dim> dilate_path_;
  vec_Vecf<Dim> center_path_;
  vec_Vecf<Dim> obs_;

  vec_E<Ellipsoid<Dim>> ellipsoids_;
  Polyhedra<Dim> polyhedra_;
  Polyhedra<Dim> intersect_polyhedra_;
  std::vector<std::shared_ptr<LineSegment>> lines_;


  bool verbose_{false};

  bool has_bounding_box_{false};
  Vecf<Dim> local_bbox_{Vecf<Dim>::Zero()};
  Vecf<Dim> global_bbox_min_{Vecf<Dim>::Zero()}; // bounding box params
  Vecf<Dim> global_bbox_max_{Vecf<Dim>::Zero()};

};
#endif
