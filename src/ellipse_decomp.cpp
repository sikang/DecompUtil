#include <decomp_util/ellipsoid_decomp.h>

EllipseDecomp::EllipseDecomp(bool verbose) {


EllipseDecomp::EllipseDecomp(const Vec3f &origin, const Vec3f &dim, bool verbose){
  has_bounding_box_ = true;
}

void EllipseDecomp::clean() {
  clear();
  obs_.clear();
}

void EllipseDecomp::clear() {
  lines_.clear();
  ellipsoids_.clear();
  polyhedrons_.clear();
  intersect_polyhedrons_.clear();
  dilate_path_.clear();
  center_path_.clear();
}

vec_LinearConstraint3f EllipseDecomp::get_constraints(){


void EllipseDecomp::add_bounding(Polyhedron &Vs) {


bool EllipseDecomp::decomp(const vec_Vec3f &poses, double offset_x) {
  clear();


void EllipseDecomp::shrink(const vec_Vec3f& path, double shrink_distance) {
  if(shrink_distance <= 0)
    return;
  for(unsigned int i = 0; i < lines_.size(); i++) {
    lines_[i]->shrink(path[i], path[i+1], shrink_distance);
    polyhedrons_[i] = lines_[i]->polyhedron();
    if(has_bounding_box_)
      add_bounding(polyhedrons_[i]);
  }
}

vec_Vec3f EllipseDecomp::cal_centers(const Polyhedra &intersect_vs) {
  vec_Vec3f path;
  for (unsigned int i = 0; i < intersect_vs.size(); i++) {
    Polyhedron vs = intersect_vs[i];

    Vec3f pt = dilate_path_[i+1];
    if(!inside_polytope(pt, vs))
      continue;

    vec_E<vec_Vec3f> bs = cal_extreme_points(vs);
    Vec3f avg = Vec3f::Zero();
    int cnt = 0;
    for (unsigned int j = 0; j < bs.size(); j++) {
      for (const auto &it : bs[j]) {
        avg += it;
        cnt++;
      }
    }
    //*** use average point as initial guess
    avg /= cnt;
#if USE_CENTROID
    Vec3f C = cal_centroid_3d(bs, avg);
#else
    Mat3f init_E;
    init_E << 0.01, 0, 0,
           0, 0.01, 0,
           0, 0, 0.01;
    Vec3f init_d = avg;
    ellipse_Vec3f e = std::make_pair(init_E, init_d);

    Vec3f C = avg;
    LinearConstraint3f Cs = cal_Axb(e.second, vs);
    if(max_ellipsoid(Cs, e)){
      C = e.second;
      ellipses_.push_back(e);
    }
#endif

    if (!inside_polytope(C, vs, 0)) // TODO: SOME WEIRD THING Happen
    {
      if(verbose_)
        std::cout << " C: " << C.transpose() << std::endl;
      path.push_back(avg);
    } else
      path.push_back(C);
  }

  /*
  double z = path.front()(2);
  for(auto &it: path)
    it(2) = z;
    */

  return path;
}

decimal_t EllipseDecomp::get_corridor_volume() {
  decimal_t V = 0;
  for (const auto &it: lines_)
    V += it->polyhedron_volume();

  for (unsigned int i = 0; i < intersect_polyhedrons_.size(); i++) {
    vec_E<vec_Vec3f> bs = cal_extreme_points(intersect_polyhedrons_[i]);
    V -= cal_volume(bs, center_path_[i+1]);
  }

  return V;
}

decimal_t EllipseDecomp::get_ellipsoid_volume() {
  decimal_t V = 0;
  for (const auto &it: lines_)
    V += it->ellipsoid_volume();
  return V;
}

