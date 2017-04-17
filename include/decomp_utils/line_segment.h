#ifndef LINE_SEGMENT_H
#define LINE_SEGMENT_H

#include <decomp_utils/data_type.h>
#include "geometry_utils.h"

class LineSegment{
  public:
    LineSegment();
    LineSegment(const Vec3f &p1, const Vec3f &p2, bool debug = false);
    void set_virtual_dim(decimal_t h, decimal_t v);
    void set_obstacles(const vec_Vec3f &obs);

    vec_Vec3f obs() const { return obs_; }
    Ellipsoid ellipsoid() const { return ellipsoid_; }
    Polyhedron polyhedron() const { return polyhedron_; }
    vec_Vec3f pts() const { return pts_; }
    decimal_t polyhedron_volume();
    decimal_t ellipsoid_volume();

    void dilate(decimal_t thr);
    void shrink(const Vec3f& p1, const Vec3f& p2, decimal_t thr);
    void adjust(pair_Vec3f& v);

  protected:
    void add_virtual_wall(Polyhedron &Vs);
    void cal_spheroid(const Vec3f& pw, const Quatf& qi, const Vec3f& center,
        Vec3f& axes, Quatf& qf);

    vec_Vec3f ps_in_polytope(const Polyhedron &Vs, const vec_Vec3f &O);
    pair_Vec3f closest_obstacle(const Ellipsoid &C, const vec_Vec3f &O);
    Ellipsoid find_ellipsoid(const Vec3f &p1, const Vec3f &p2);
    Polyhedron find_polyhedron(const Ellipsoid &E);

    // Input:
    Vec3f p1_;
    Vec3f p2_;
    vec_Vec3f obs_;

    // Output:
    Ellipsoid ellipsoid_;
    Polyhedron polyhedron_;
    vec_Vec3f pts_;

    // Debug Mode:
    bool debug_;

    Vec3f min_; // bounding box for visualization
    Vec3f max_;
    Vec3f origin_;
    Vec3f dim_;

    decimal_t robot_radius_;
    decimal_t VIRTUAL_H = 5;
    decimal_t VIRTUAL_V = 2;
};
#endif
