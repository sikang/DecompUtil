#ifndef ELLIPSE_DECOMP_DATA_TYPE_H
#define ELLIPSE_DECOMP_DATA_TYPE_H

#include <data_utils.h>

//typedef std::pair<Vec3f, Vec3f> pair_Vec3f;
typedef vec_E<vec_E<vec_Vec3f>> bound_Vec3f; // compose by extreme points
typedef std::pair<Mat3f, Vec3f> Ellipsoid;
typedef vec_E<Ellipsoid> vec_Ellipsoid;

constexpr decimal_t epsilon_ = 1e-6; // numerical calculation effect
#endif
