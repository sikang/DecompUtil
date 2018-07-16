/**
 * @file line_segment.h
 * @brief LineSegment Class
 */
#ifndef LINE_SEGMENT_H
#define LINE_SEGMENT_H

#include <decomp_util/data_type.h>
#include <decomp_util/geometry_utils.h>

/**
 * @brief Line Segment Class
 *
 * The basic element in EllipsoidDecomp
 */
template <int Dim>
class LineSegment {
  public:
    ///Simple constructor
    LineSegment();
    /**
     * @brief Basic constructor
     * @param p1 One end of the line seg
     * @param p2 The other end of the line seg
     */
    LineSegment(const Vecf<Dim> &p1, const Vecf<Dim> &p2);
    /**
     * @brief Adding local bounding box around line seg
     * @param Dim Distance in corresponding axis
     *
     * This virtual bounding box is parallel to the line segment, the x,y,z axes are not w.r.t the world coordinate system, but instead, x-axis is parallel to the line, y-axis is perpendicular to the line and world z-axis, z-axis is perpendiculat to the line and y-axis
     */
    void set_local_bbox(const Vecf<Dim>& bbox);
    ///Import obstacle points
    void set_obstacles(const vec_Vecf<Dim> &obs);
    ///Get obstacel points
    vec_Vecf<Dim> get_obs() const { return obs_; }
    ///Get ellipsoid
    Ellipsoid<Dim> get_ellipsoid() const { return ellipsoid_; }
    ///Retrieve polyhedron
    Polyhedron<Dim> get_polyhedron() const { return polyhedron_; }
    ///Calculate the volume of polyhedron
    decimal_t polyhedron_volume();
    ///Calculate the volume of ellipsoid
    decimal_t ellipsoid_volume();

    /**
     * @brief Infalte the line segment
     * @param radius the offset added to the long semi-axis
     */
    void dilate(decimal_t radius = 0);
    /**
     * @brief Shrink the polyhedron
     * @param p1 One end of the line seg
     * @param p2 The other end of the line seg
     * @param shrink_distance Shrink distance
     */
    void shrink(const Vecf<Dim>& p1, const Vecf<Dim>& p2, double shrink_distance);
    /**
     * @brief Adjust the norm of half plane to prevent cutting off the line seg whhile shrinking
     *
     * Details are introduced in the related paper
     */
    void adjust(Face& v);

  protected:
    void add_local_bbox(Polyhedron<Dim> &Vs);
    void cal_spheroid(const Vec3f& pw, const Quatf& qi, const Vec3f& center,
        Vec3f& axes, Quatf& qf);

    Ellipsoid<Dim> find_ellipsoid(const Vecf<Dim> &p1, const Vecf<Dim> &p2, double offset_x = 0);
    Polyhedron<Dim> find_polyhedron(const Ellipsoid<Dim> &E);

    /// One end of line segment, input
    Vecf<Dim> p1_;
    /// The other end of line segment, input
    Vecf<Dim> p2_;
    /// Obstacles, input
    vec_Vecf<Dim> obs_;

    /// Output ellipsoid
    Ellipsoid<Dim> ellipsoid_;
    /// Output polyhedron
    Polyhedron<Dim> polyhedron_;

    /// Local bounding box along the line segment
    Vecf<Dim> local_bbox_{Vecf<Dim>::Zero()};
};
#endif
