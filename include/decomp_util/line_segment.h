/**
 * @file line_segment.h
 * @brief LineSegment Class
 */
#ifndef LINE_SEGMENT_H
#define LINE_SEGMENT_H

#include <decomp_util/decomp_base.h>

/**
 * @brief Line Segment Class
 *
 * The basic element in EllipsoidDecomp
 */
template <int Dim>
class LineSegment : public DecompBase<Dim> {
  public:
    ///Simple constructor
    LineSegment() {};
    /**
     * @brief Basic constructor
     * @param p1 One end of the line seg
     * @param p2 The other end of the line seg
     */
    LineSegment(const Vecf<Dim> &p1, const Vecf<Dim> &p2) : p1_(p1), p2_(p2) {}
    /**
     * @brief Infalte the line segment
     * @param radius the offset added to the long semi-axis
     */
    void dilate(decimal_t radius = 0) {
      find_ellipsoid(radius);
      this->find_polyhedron();
      add_local_bbox(this->polyhedron_);
    }
    /**
     * @brief Adjust the norm of half plane to prevent cutting off the line seg while shrinking
     *
     * Details are introduced in the related paper
     */
    void adjust(Face& v);

  protected:
    ///Add the bounding box
    void add_local_bbox(Polyhedron<Dim> &Vs) {
      if(this->local_bbox_.norm() == 0)
        return;
      //**** virtual walls parallel to path p1->p2
      Vecf<Dim> dir = (p2_ - p1_).normalized();
      Vecf<Dim> dir_h = Vecf<Dim>::Zero();
      dir_h(0) = dir(1), dir_h(2) = -dir(0);
      if (dir_h.norm() == 0) {
        if(Dim == 2)
          dir_h << -1, 0;
        else
          dir_h << -1, 0, 0;
      }
      dir_h = dir_h.normalized();

      // along x
      Vecf<Dim> pp1 = p1_ + dir_h * local_bbox_(1);
      Vecf<Dim> pp2 = p1_ - dir_h * local_bbox_(1);
      Vs.add(Hyperplane<Dim>(pp1, dir_h));
      Vs.add(Hyperplane<Dim>(pp2, -dir_h));

      // along y
      Vecf<Dim> pp3 = p2_ + dir * local_bbox_(0);
      Vecf<Dim> pp4 = p1_ - dir * local_bbox_(0);
      Vs.add(Hyperplane<Dim>(pp3, dir));
      Vs.add(Hyperplane<Dim>(pp4, -dir));

      // along z
      if(Dim > 2) {
        Vecf<Dim> dir_v = dir.cross(dir_h);
        Vecf<Dim> pp5 = p1_ + dir_v * local_bbox_(2);
        Vecf<Dim> pp6 = p1_ - dir_v * local_bbox_(2);
        Vs.add(Hyperplane<Dim>(pp5, dir_v));
        Vs.add(Hyperplane<Dim>(pp6, -dir_v));
      }
    }

    void find_ellipsoid(double offset_x) {
      const decimal_t f = (p1_ - p2_).norm() / 2;
      Matf<Dim, Dim> C = f * Matf<Dim, Dim>::Identity();
      Vecf<Dim> axes = f * Vecf<Dim>::One();
      C(0, 0) += offset_x;
      axes(0) += offset_x;

      if(axes(0) > 0) {
        double ratio = axes(1) / axes(0);
        axes *= ratio;
        C *= ratio;
      }

      const auto Ri = vec_to_rotation(p2_ - p1_);
      C = Ri * C * Ri.transpose();

      Ellipsoid<Dim> E(C, (p1_ + p2_) / 2);
      auto Rf = Ri;

      auto obs = E.points_inside(obs_);
      auto obs_inside = obs;
      //**** decide short axes
      while (!obs_inside.empty()) {
        const auto pw = E.closest_point(obs_inside);
        auto p = Ri.transpose() * (pw - E.d()); // to ellipsoid frame
        if(Dim > 2) {
          const decimal_t roll = atan2(p(2), p(1));
          Rf = Ri * Quatf(cos(roll / 2), sin(roll / 2), 0, 0);
          p = Rf.inverse() * (pw - E.second);
        }

        axes(1) = std::abs(p(1)) / std::sqrt(1 - std::pow(p(0) / axes(0), 2));
        auto new_C = Matf<Dim, Dim>::Identity();
        new_C(0, 0) = axes(0);
        new_C(1, 1) = axes(1);
        if(Dim > 2)
          new_C(2, 2) = axes(1);
        E.C_ = Rf * new_C * Rf.transpose();

        vec_Vecf<Dim> obs_new;
        for(const auto &it: obs_inside) {
          if(E.dist(it) != 1)
            obs_new.push_back(it);
        }
        obs_inside = obs_new;
      }

      if(Dim == 2) {
        this->ellipsoid_ = E;
        return;
      }

      //**** reset ellipsoid with old axes(2)
      C = f * Matf<Dim, Dim>::Identity();
      C(0, 0) = axes(0);
      C(1, 1) = axes(1);
      C(2, 2) = axes(2);
      E.C_ = Rf * C * Rf.transpose();
      obs_inside = E.points_inside(obs);

      while (!obs_inside.empty()) {
        const auto pw = E.closest_point(obs_inside);

        Vec3f p = Rf.transpose() * (pw - E.d());
        axes(2) =
          std::abs(p(2)) / sqrt(1 - pow(p(0) / axes(0), 2) - pow(p(1) / axes(1), 2));
        auto new_C = Matf<Dim, Dim>::Identity();
        new_C(0, 0) = axes(0);
        new_C(1, 1) = axes(1);
        new_C(2, 2) = axes(2);
        E.C_ = Rf * new_C * Rf.transpose();

        vec_Vecf<Dim> obs_new;
        for(const auto &it: obs_inside) {
          if(E.dist(it) != 1)
            obs_new.push_back(it);
        }
        obs_inside = obs_new;
      }

      this->ellipsoid_ = E;
    }

    /// One end of line segment, input
    Vecf<Dim> p1_;
    /// The other end of line segment, input
    Vecf<Dim> p2_;
};
#endif
