#include <decomp_util/seed_decomp.h>

SeedDecomp::SeedDecomp(const Vec3f& p): p_(p) {}

void SeedDecomp::set_virtual_dim(decimal_t x, decimal_t y, decimal_t z){
  virtual_x_ = x;
  virtual_y_ = y;
  virtual_z_ = z;
}

void SeedDecomp::set_obstacles(const vec_Vec3f& obs) {
  Polyhedron vs;
  add_virtual_wall(vs);
  obs_ = ps_in_polytope(vs, obs);
}

void SeedDecomp::dilate(decimal_t radius){
  dilate(Vec3f(radius, radius, radius), Mat3f::Identity());
}

void SeedDecomp::dilate(const Vec3f& axes, const Mat3f& R){
  Mat3f C = Mat3f::Identity();
  C(0, 0) = axes(0);
  C(1, 1) = axes(1);
  C(2, 2) = axes(2);
  C = R * C * R.transpose();
  Ellipsoid E = std::make_pair(C, p_);
  dilate(E);
}

void SeedDecomp::dilate(const Ellipsoid& E) {
  double radius = 1.0; // default radius is 1
  if(!obs_.empty()) {
    Face v = closest_obstacle(E, obs_);
    decimal_t b = v.p.dot(v.n);
    radius = std::abs(b - v.n.dot(p_)) / std::sqrt((v.n.dot(E.first * E.first.transpose()*v.n)));
    //printf("radius = %f\n", radius);
  }

  ellipsoid_.first = radius * E.first;
  ellipsoid_.second = E.second;

  //**** find half-space
  Polyhedron Vs;
  vec_Vec3f O_remain = obs_;
  while (!O_remain.empty()) {
    Face v = closest_obstacle(ellipsoid_, O_remain);
    Vs.push_back(v);
    Vec3f a = v.n;
    decimal_t b = v.p.dot(a);
    vec_Vec3f O_tmp;
    for (const auto &it : O_remain) {
      if (a.dot(it) - b < 0)
        O_tmp.push_back(it);
    }
    O_remain = O_tmp;
  }

  add_virtual_wall(Vs);
  polyhedron_ = Vs;
}


void SeedDecomp::shrink(decimal_t thr) {
}


