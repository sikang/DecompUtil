#include <decomp_util/ellipse_decomp.h>

EllipseDecomp::EllipseDecomp(const Vec3f &origin, const Vec3f &dim, bool verbose){
  min_ = origin;
  max_ = origin + dim;
  verbose_ = verbose;
  if(verbose_){
    printf(ANSI_COLOR_GREEN "DECOMP VERBOSE ON! \n");
    info();
  }
}

void EllipseDecomp::info() {
    printf("Min: [%f, %f, %f]\n", min_(0), min_(1), min_(2));
    printf("Max: [%f, %f, %f]\n" ANSI_COLOR_RESET, max_(0), max_(1), max_(2));
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
  vec_LinearConstraint3f constraints;
  constraints.resize(polyhedrons_.size());
  for (unsigned int i = 0; i < polyhedrons_.size(); i++){
    Vec3f pt = (center_path_[i] + center_path_[i+1])/2;
    constraints[i] = cal_Axb(pt, polyhedrons_[i]);
  }
  return constraints;
}

Polyhedra EllipseDecomp::get_polyhedra(int id) const {
  Polyhedra polys;
  if(id == 0)
    return polyhedrons_;
  else
    return intersect_polyhedrons_;
}

vec_Vec3f EllipseDecomp::get_pts() const {
  vec_Vec3f pts;
  for(const auto& l: lines_) {
    vec_Vec3f ps = l->pts();
    pts.insert(pts.end(), ps.begin(), ps.end());
  }

  return pts;
}

void EllipseDecomp::add_bounding(Polyhedron &Vs) {
  //**** add bound along X, Y, Z axis
  //*** Z
  Vs.push_back(Face(Vec3f(0, 0, max_(2)), Vec3f(0, 0, 1)));
  Vs.push_back(Face(Vec3f(0, 0, min_(2)), Vec3f(0, 0, -1)));

  //*** X
  Vs.push_back(Face(Vec3f(max_(0), 0, 0), Vec3f(1, 0, 0)));
  Vs.push_back(Face(Vec3f(min_(0), 0, 0), Vec3f(-1, 0, 0)));
  //*** Y
  Vs.push_back(Face(Vec3f(0, max_(1), 0), Vec3f(0, 1, 0)));
  Vs.push_back(Face(Vec3f(0, min_(1), 0), Vec3f(0, -1, 0)));
}

bool EllipseDecomp::decomp(const vec_Vec3f &poses,
                           decimal_t d,
                           decimal_t ds,
                           decimal_t h) {
  clear();
  if (poses.size() < 2)
  {
    if(verbose_)
      printf(ANSI_COLOR_RED "Decomp failed, poses size: %zu\n" ANSI_COLOR_RESET, poses.size());
    return false;
  }

  unsigned int N = poses.size() - 1;
  lines_.resize(N);
  ellipsoids_.resize(N);
  polyhedrons_.resize(N);
  intersect_polyhedrons_.resize(N-1);

  for (unsigned int i = 0; i < N; i++) {
    lines_[i] = std::make_shared<LineSegment>(poses[i], poses[i+1], false);
    lines_[i]->set_virtual_dim(10, h);
    lines_[i]->set_obstacles(obs_);
    lines_[i]->dilate(d);

    ellipsoids_[i] = lines_[i]->ellipsoid();
    polyhedrons_[i] = lines_[i]->polyhedron();
  }

  for (unsigned int i = 0; i < poses.size() - 2; i++){
    intersect_polyhedrons_[i] = polytope_intersection(polyhedrons_[i], polyhedrons_[i+1]);
    add_bounding(intersect_polyhedrons_[i]);
  }

  dilate_path_ = poses;
  center_path_ = cal_centers(intersect_polyhedrons_);
  center_path_.push_back(poses.back());
  center_path_.insert(center_path_.begin(), poses.front());

  for(unsigned int i = 0; i < lines_.size(); i++){
    //lines_[i]->shrink(dilate_path_[i], dilate_path_[i+1], d);
    lines_[i]->shrink(center_path_[i], center_path_[i+1], ds * d);
    polyhedrons_[i] = lines_[i]->polyhedron();
    add_bounding(polyhedrons_[i]);
  }
  return true;
}


//**** prev_v is the second last point
void EllipseDecomp::change_end(decimal_t z_up, decimal_t z_down){
  const Vec3f p1 = center_path_[center_path_.size() - 2];
  const Vec3f p2 = center_path_.back();
  const Vec3f n = (p2 - p1).normalized();
  Face p(p2, n);
  Polyhedron vs = polyhedrons_.back();

  //**** virtual walls parallel to path p1->p2
  Vec3f dir = p2 - p1;
  dir /= dir.norm();
  Vec3f dir_h(-dir(1), dir(0), 0);
  if (dir_h == Vec3f::Zero())
    dir_h << -1, 0, 0;
  Vec3f pp1 = p1 + dir_h * 3.0;
  Vec3f pp2 = p1 - dir_h * 3.0;
  vs.push_back(Face(pp1, dir_h));
  vs.push_back(Face(pp2, -dir_h));

  Vec3f dir_v = dir.cross(dir_h);
  Vec3f pp3 = p1 + dir_v * z_up;
  Vec3f pp4 = p1 - dir_v * z_down;
  vs.push_back(Face(pp3, dir_v));
  vs.push_back(Face(pp4, -dir_v));
  vec_Vec3f valid_pts = plane_polytope_intersection(p, vs);
  Vec3f c = cal_centroid_2d(valid_pts, p);
  if (inside_polytope(c, polyhedrons_.back())) {
    center_path_.back() = c;
    if(verbose_)
      printf(ANSI_COLOR_GREEN "Change end! \n" ANSI_COLOR_RESET);
  }
  else{
    if(verbose_)
      printf(ANSI_COLOR_RED "Change end failed! \n" ANSI_COLOR_RESET);
  }
}


vec_Vec3f EllipseDecomp::cal_centers(const Polyhedra &intersect_vs) {
  vec_Vec3f path;
  for (unsigned int i = 0; i < intersect_vs.size(); i++) {
    Polyhedron vs = intersect_vs[i];

    Vec3f pt = dilate_path_[i+1];
    if(!inside_polytope(pt, vs))
      continue;

    Vec3f dir_v(0, 0, 1);
    Vec3f pp3 = pt + dir_v * 0.2;
    Vec3f pp4 = pt - dir_v * 0.2;
    vs.push_back(Face(pp3, dir_v));
    vs.push_back(Face(pp4, -dir_v));

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

