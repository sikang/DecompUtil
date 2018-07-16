/**
 * @file geometric_utils.h
 * @brief basic geometry utils
 */
#ifndef DECOMP_GEOMETRIC_UTILS_H
#define DECOMP_GEOMETRIC_UTILS_H

#include <iostream>
#include <decomp_util/data_utils.h>
#include <Eigen/Eigenvalues>

///Compensate for numerical error
constexpr decimal_t epsilon_ = 1e-6; // numerical calculation error

/// Calculate eigen values
template <int Dim>
Vecf<Dim> eigen_value(const Matf<Dim, Dim>& A) {
  Eigen::SelfAdjointEigenSolver<Matf<Dim, Dim>> es(A);
  return es.eigenvalues();
}

/// Sort points on the same plane in the counter-clockwise order
vec_Vecf<Dim> sort_pts(const vec_Vecf<Dim> &pts) {
  /// if empty, dont sort
  if(pts.empty())
    return pts;
  /// calculate center point
  Vecf<Dim> avg = Vecf<Dim>::Zero();
  for (const auto& pt : pts)
    avg += pt;
  avg /= pts.size();

  /// sort in body frame
  vec_E<std::pair<decimal_t, Vecf<Dim>>> pts_valued;
  pts_valued.resize(pts.size());
  for (unsigned int i = 0; i < pts.size(); i++) {
    decimal_t theta = atan2(pts[i](1) - avg(1), pts[i](0) - avg(0));
    pts_valued[i] = std::make_pair(theta, pts[i]);
  }

  std::sort(pts_valued.begin(), pts_valued.end(),
            [](const std::pair<decimal_t, Vecf<Dim>> &i,
               const std::pair<decimal_t, Vecf<Dim>> &j) {
              return i.first < j.second;});
  vec_Vecf<Dim> pts_sorted(pts_valued.size());
  for (size_t i = 0; i < pts_valued.size(); i++)
    pts_sorted[i] = pts_valued[i].second;
  return pts_sorted;
}





///Find the closest point
template <int Dim>
Vecf<Dim> closest_point(const Ellipsoid &E,
                        const vec_Vec3f &O, int& id) {
  int cnt = 0;
  Vecf<Dim> pt = Vecf<Dim>::Zero();
  decimal_t dist = std::numeric_limits<decimal_t>::max();
  for (const auto &it : O) {
    decimal_t d = (E.first.inverse() * (it - E.second)).norm();
    if (d < dist) {
      dist = d;
      pt = it;
      id = cnt;
    }
    cnt ++;
  }
  return pt;
}

///Find the closest hyperplane from the closest point
template <int Dim>
Hyperplane<Dim> closest_hyperplance(const Ellipsoid<Dim> &E,
                                    const vec_Vecf<Dim> &O) {
  int id = -1;
  const auto closest_pt = closest_point(E, O, id);
  if(id >= 0) {
    const auto n = E.first.inverse() * E.first.inverse().transpose() *
      (closest_pt - E.second);
    return std::make_pair(closest_pt, n.normalized());
  }
  else
    return Hyperplane<Dim>();
}

/// Calculate points inside ellipsoid
template <int Dim>
vec_Vecf<Dim> points_inside_ellipsoid(const Ellipsoid<Dim> &E,
                                      const vec_Vecf<Dim> &O) {
  vec_Vecf<Dim> new_O;
  for (const auto &it : O) {
    decimal_t d = (E.first.inverse() * (it - E.second)).norm();
    if (d < 1)
      new_O.push_back(it);
  }
  return new_O;
}

//**** Find intersection between two Line
//*** return false if they are not intersected
bool intersect(const pair_Vec3f &v1, const pair_Vec3f &v2, Vec3f &pi);

//**** Find intersection between multiple lines
vec_Vec3f line_intersects(const vec_E<pair_Vec3f> &lines);

//**** Find extreme points of polytope
vec_E<vec_Vec3f> cal_extreme_points(const Polyhedron &vts);

//**** Find intersect polygon between a plane and polytope
vec_Vec3f plane_polytope_intersection(const Face &plane,
                                      const Polyhedron& vts);

//**** uniformly sample path into many segments
vec_Vec3f path_downsample(const vec_Vec3f& ps, decimal_t d);
vec_Vec3f path_downsample_i(const vec_Vec3f& ps, int cnt);

//**** crop path
vec_Vec3f path_crop(const vec_Vec3f& path, decimal_t d);

//**** find the intersection of two polytopes
Polyhedron polytope_intersection(const Polyhedron& vs1,
    const Polyhedron& vs2);

//**** Create triangles from a face
vec_E<vec_Vec3f> chop_triangle(const vec_Vec3f& pts);

//**** Calculate the volume of a polytope
decimal_t cal_volume(const vec_E<vec_Vec3f>& fs, const Vec3f& pt_inside);

//**** Find the centroid of a polygon
Vec3f cal_centroid_2d(const vec_Vec3f &pts, const Face&p);


//**** Calculate centroid of a polytope
Vec3f cal_centroid_3d(const vec_E<vec_Vec3f>& fs, const Vec3f& pt_inside);

decimal_t cal_closest_dist(const Vec3f& pt, const Polyhedron& vs);


Quatf vec_to_quaternion(const Vec3f &v);

bool inside_ellipsoid(const Ellipsoid &E,
                      const vec_Vec3f &O);

bool closest_pt(const Ellipsoid &E,
    const vec_Vec3f &O,
    Vec3f &best_v, int& id);

vec_Vec3f ps_in_ellipsoid(const Ellipsoid &E,
                          const vec_Vec3f &O);

#endif
