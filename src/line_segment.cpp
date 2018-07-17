#include <decomp_util/line_segment.h>

template <int Dim>
Ellipsoid<Dim> LineSegment<Dim>::find_ellipsoid(double offset_x) {
template <int Dim>
void LineSegment<Dim>::add_local_bbox(Polyhedron<Dim> &Vs) {
}

template <int Dim>
void LineSegment<Dim>::dilate(decimal_t radius) {


template <int Dim>
void LineSegment<Dim>::shrink(const Vec3f& p1, const Vec3f& p2, double shrink_distance) {
  /*
  for (auto &it : polyhedron_) {
    decimal_t b = it.p.dot(it.n);
    decimal_t d1 = it.n.dot(p1) - b;
    decimal_t d2 = it.n.dot(p2) - b;
    decimal_t d = -std::max(d1, d2) - 0.1;
    d = d < shrink_distance ? d : shrink_distance;
    if (d > 0.01)
      it.p -= d * it.n;
  }
  */
}
