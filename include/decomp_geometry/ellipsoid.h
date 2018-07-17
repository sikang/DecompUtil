/**
 * @file ellipsoid.h
 * @brief Ellipsoid class
 */

#ifndef DECOMP_ELLIPSOID_H
#define DECOMP_ELLIPSOID_H

#include <decomp_basis/data_type.h>
#include <decomp_geometry/polyhedron.h>

template <int Dim>
struct Ellipsoid {
  Ellipsoid() {}
  Ellipsoid(const Matf<Dim, Dim>& C, const Vecf<Dim>& d) : C_(C), d_(d) {}

  /// Calculate distance to the center
  decimal_t dist(const Vecf<Dim>& pt) {
    return (C_.inverse() * (pt - d_)).norm();
  }

  /// Check if the point is inside, non-exclusive
  bool inside(const Vecf<Dim>& pt) {
      return dist(pt) <= 1;
  }

  /// Calculate points inside ellipsoid, non-exclusive
  vec_Vecf<Dim> points_inside(const vec_Vecf<Dim> &O) {
    vec_Vecf<Dim> new_O;
    for (const auto &it : O) {
      if (inside(it))
        new_O.push_back(it);
    }
    return new_O;
  }

  ///Find the closest point
  Vecf<Dim> closest_point(const vec_Vecf<Dim> &O) {
    Vecf<Dim> pt = Vecf<Dim>::Zero();
    decimal_t min_dist = std::numeric_limits<decimal_t>::max();
    for (const auto &it : O) {
      decimal_t d = dist(it);
      if (d < min_dist) {
        min_dist = d;
        pt = it;
      }
    }
    return pt;
  }

  ///Find the closest hyperplane from the closest point
  Hyperplane<Dim> closest_hyperplane(const vec_Vecf<Dim> &O) {
    const auto closest_pt = closest_point(O);
    const auto n = C_.inverse() * C_.inverse().transpose() *
      (closest_pt - d_);
    return Hyperplane<Dim>(closest_pt, n.normalized());
  }

  vec_Vecf<Dim> sample(int num) const {
    vec_Vecf<Dim> pts;
    if(Dim == 2) {
      decimal_t dyaw = M_PI*2/num;
      for(decimal_t yaw = 0; yaw < M_PI*2; yaw+=dyaw) {
        Vecf<Dim> pt;
        pt << cos(yaw), sin(yaw);
        pts.push_back(C_ * pt + d_);
      }
    }
    return pts;
  }

  /// Get ellipsoid volume
  decimal_t volume() const {
    return C_.determinant();
  }

  /// Get C matrix
  Matf<Dim, Dim> C() const {
    return C_;
  }

  /// Get center
  Vecf<Dim> d() const {
    return d_;
  }

  Matf<Dim, Dim> C_;
  Vecf<Dim> d_;
};

typedef Ellipsoid<2> Ellipsoid2D;

typedef Ellipsoid<3> Ellipsoid3D;

#endif
