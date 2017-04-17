#ifndef ITERATIVE_DECOMP_H
#define ITERATIVE_DECOMP_H

#include <decomp_utils/ellipse_decomp.h>

class IterativeDecomp : public EllipseDecomp
{
  public:
    IterativeDecomp(const Vec3f &origin, const Vec3f &dim, bool verbose = false);
    bool decomp_iter(const vec_Vec3f &poses,
                     decimal_t d = 0.0,
                     decimal_t ds = 1.0,
                     decimal_t h = 3.0,
                     int iter_num = 5,
		     bool fixed = false);
  private:
    vec_Vec3f simplify(const vec_Vec3f& path);

};
#endif
