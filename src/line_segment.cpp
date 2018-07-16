#include <decomp_util/line_segment.h>

template <int Dim>
LineSegment<Dim>::LineSegment(const Vecf<Dim> &p1, const Vecf<Dim> &p2):
  p1_(p1), p2_(p2) {}

template <int Dim>
void LineSegment<Dim>::set_local_bbox(const Vecf<Dim>& bbox) {
  local_bbox_ = bbox;
}

template <int Dim>
void LineSegment<Dim>::set_obstacles(const vec_Vecf<Dim>& obs) {
  // only consider points inside local bbox
  Polyhedron<Dim> vs;
  add_local_bbox(vs);
  obs_ = points_inside_polyhedron(vs, obs);
}

template <int Dim>
Ellipsoid<Dim> LineSegment<Dim>::find_ellipsoid(const Vecf<Dim> &p1,
                                                const Vecf<Dim> &p2,
                                                double offset_x) {
  const decimal_t f = (p1 - p2).norm() / 2;
  Matf<Dim, Dim> C = f * Matf<Dim, Dim>::Identity();
  Vecf<Dim> axes = f * Vecf<Dim>::One();
  C(0, 0) += offset_x;
  axes(0) += offset_x;

  if(axes(0) > 0) {
    double ratio = axes(1) / axes(0);
    axes *= ratio;
    C *= ratio;
  }

  const auto Ri = vec_to_rotation(p2 - p1);
  C = Ri * C * Ri.transpose();

  Ellipsoid<Dim> E = std::make_pair(C, (p1 + p2) / 2);
  auto Rf = Ri;

  vec_Vecf<Dim> obs = points_inside_ellipsoid(E, obs_);
  Vecf<Dim> prev_pw;
  bool removed = false;
  //**** decide short axes
  while (inside_ellipsoid(E, obs)) {
    int id = -1;
    const auto pw = closest_point(E, obs, id);
    auto p = Ri.transpose() * (pw - E.second); // to ellipsoid frame
    if(Dim > 2) {
      const decimal_t roll = atan2(p(2), p(1));
      Rf = Ri * Quatf(cos(roll / 2), sin(roll / 2), 0, 0);
      p = Rf.inverse() * (pw - E.second);
    }

    axes(1) = std::abs(p(1)) / std::sqrt(1 - std::pow(p(0) / axes(0), 2));
    E.first = Matf<Dim, Dim>::Identity();
    E.first(0, 0) = axes(0);
    E.first(1, 1) = axes(1);
    if(Dim > 2)
      E.first(2, 2) = axes(1);
    E.first = Rf * E.first * Rf.transpose();

    obs.erase(obs.begin() + id); // remove pw
    if(removed)
      obs.push_back(prev_pw);
    removed = true;
    prev_pw = pw;
  }

  if(Dim == 2)
    return E;

  //**** reset ellipsoid with old axes(2)
  E.first = f * Matf<Dim, Dim>::Identity();
  E.first(0, 0) = axes(0);
  E.first(1, 1) = axes(1);
  E.first(2, 2) = axes(2);
  E.first = Rf * E.first * Rf.transpose();
  obs = points_inside_ellipsoid(E, obs);

  removed = false;
  while (inside_ellipsoid(E, obs)) {
    Vecf<Dim> pw = Vecf<Dim>::Zero();
    int id = -1;
    closest_pt(E, obs, pw, id);

    Vec3f p = qf.inverse() * (pw - E.second);
    axes(2) =
      fabs(p(2)) / sqrt(1 - pow(p(0) / axes(0), 2) - pow(p(1) / axes(1), 2));
    E.first = Mat3f::Identity();
    E.first(0, 0) = axes(0);
    E.first(1, 1) = axes(1);
    E.first(2, 2) = axes(2);
    E.first = qf * E.first * qf.conjugate();

    Os.erase(Os.begin() + id); // remove pw
    if(removed)
      Os.push_back(prev_pw);
    removed = true;
    prev_pw = pw;
  }

  return E;
}

Polyhedron LineSegment::find_polyhedron(const Ellipsoid& E){
  //**** find half-space
  Polyhedron Vs;
  vec_Vec3f O_remain = obs_;
  while (!O_remain.empty()) {
    Face v = closest_obstacle(E, O_remain);
    //adjust(v);
    Vs.push_back(v);
    Vec3f a = v.n;
    decimal_t b = v.p.dot(a);
    vec_Vec3f O_tmp;
    for (const auto &it : O_remain) {
      if (a.dot(it) - b < 0)
        O_tmp.push_back(it);
    }
    O_remain = O_tmp;
    /*
    std::cout << "a: " << a.transpose() << std::endl;
    std::cout << "b: " << b << std::endl;
    */
  }

  return Vs;
}

template <int Dim>
void LineSegment<Dim>::add_local_bbox(Polyhedron<Dim> &Vs) {
  if(local_bbox_.norm() == 0)
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
  Vs.push_back(std::make_pair(pp1, dir_h));
  Vs.push_back(std::make_pair(pp2, -dir_h));

  // along y
  Vecf<Dim> pp3 = p2_ + dir * local_bbox_(0);
  Vecf<Dim> pp4 = p1_ - dir * local_bbox_(0);
  Vs.push_back(std::make_pair(pp3, dir));
  Vs.push_back(std::make_pair(pp4, -dir));

  if(Dim > 2) {
    Vecf<Dim> dir_v = dir.cross(dir_h);
    Vecf<Dim> pp5 = p1_ + dir_v * local_bbox_(2);
    Vecf<Dim> pp6 = p1_ - dir_v * local_bbox_(2);
    Vs.push_back(std::make_pair(pp5, dir_v));
    Vs.push_back(std::make_pair(pp6, -dir_v));
  }
}

void LineSegment::dilate(decimal_t radius) {
  ellipsoid_ = find_ellipsoid(p1_, p2_, radius);
  polyhedron_ = find_polyhedron(ellipsoid_);
  add_local_bbox(polyhedron_);
}


decimal_t LineSegment::ellipsoid_volume() {
  return ellipsoid_.first.determinant();
}

decimal_t LineSegment::polyhedron_volume() {
  vec_E<vec_Vec3f> bs = cal_extreme_points(polyhedron_);
  return cal_volume(bs, (p1_ + p2_)/2);
}


void LineSegment::shrink(const Vec3f& p1, const Vec3f& p2, double shrink_distance) {
  for (auto &it : polyhedron_) {
    decimal_t b = it.p.dot(it.n);
    decimal_t d1 = it.n.dot(p1) - b;
    decimal_t d2 = it.n.dot(p2) - b;
    decimal_t d = -std::max(d1, d2) - 0.1;
    d = d < shrink_distance ? d : shrink_distance;
    if (d > 0.01)
      it.p -= d * it.n;
  }
}
