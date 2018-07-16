#include <decomp_util/data_type.h>

template <int Dim>
struct Ellipsoid {
  Ellipsoid() {}
  Ellipsoid(const Matf<Dim, Dim>& C, const Vecf<Dim>& d) : C_(C), d_(d) {}

  /// Check if the point is inside, non-exclusive
  bool inside(const Vecf<Dim>& pt) {
      decimal_t d = (C_.inverse() * (it - d_)).norm();
      return d <= 1;
  }

  Matf<Dim, Dim> C_;
  Vecf<Dim> d_;
};

typedef Ellipsoid<2> Ellipsoid2D;

typedef Ellipsoid<3> Ellipsoid3D;
