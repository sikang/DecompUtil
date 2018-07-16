#include <decomp_util/data_type.h>

///Hyperplane class
template <int Dim>
struct Hyperplane {
  Hyperplane() {}
  Hyperplane(const Vecf<Dim>& p, const Vecf<Dim>& n) : p_(p), n_(n) {}

  /// Calculate the signed distance from point
  decimal_t signed_dist(const Vecf<Dim>& p) {
    return n_.dot(p - p_);
  }

  /// Calculate the distance from point
  decimal_t dist(const Vecf<Dim>& p) {
    return std::abs(signed_dist(p));
  }

  /// Point on the plane
  Vecf<Dim> p_;
  /// Normal of the plane, directional
  Vecf<Dim> n_;
};

///Hyperplane2D: first is the point on the hyperplane, second is the normal
typedef Hyperplane<2> Hyperplane2D
///Hyperplane3D: first is the point on the hyperplane, second is the normal
typedef Hyperplane<3> Hyperplane3D


///Polyhedron class
template <int Dim>
struct Polyhedron {
  Polyhedron() {}
  Polyhedron(const vec_E<Hyperplane<Dim>>& vs) : vs_(vs) {}

  /// Check if the point is inside polyhedron, non-exclusive
  bool inside(const Vecf<Dim>& pt) {
    for (const auto& v : vs_) {
      if (v.signed_dist(pt) > 0)
        return false;
    }
    return true;
  }

  /// Calculate points inside polyhedron
  vec_Vecf<Dim> points_inside(const vec_Vecf<Dim> &O) {
    vec_Vecf<Dim> new_O;
    for (const auto &it : O) {
      if (inside(it))
        new_O.push_back(it);
    }
    return new_O;
  }

  /// Calculate normals, used for visualization
  vec_E<std::pair<Vecf<Dim>, Vecf<Dim>> cal_normals() {
    vec_E<std::pair<Vecf<Dim>, Vecf<Dim>> ns(vs_.size());
    for (size_t i = 0; i < vs_.size(); i++)
      ns[i] = std::make_pair(it.p_, it.n_); // fist is point, second is normal
    return ns;
  }

  /// Hyperplane array
  vec_E<Polyhedron<Dim>> vs_; // normal must go outside

};

///Polyhedron2D, consists of 2D hyperplane
typedef Polyhedron<2> Polyhedron2D;
///Polyhedron3D, consists of 3D hyperplane
typedef Polyhedron<3> Polyhedron3D;

///Vector of Polyhedron
template <int Dim>
using Polyhedra = vec_E<Polyhedron<Dim>>;
///Vector of Polyhedron2D
typedef Polyhedra<2> Polyhedra2D;
///Vector of Polyhedron3D
typedef Polyhedra<3> Polyhedra3D;

///[A, b] for \f$Ax < b\f$
template <int Dim>
struct LinearConstraint {
  ///Null constructor
  LinearConstraint() {}
  /// Construct from \f$A, b\f$ directly, s.t \f$Ax < b\f$
  LinearConstraint(const MatDNf<Dim>& A, const VeDf& b) : A_(A), b_(b) {}
  /**
   * @brief Construct from a inside point and hyperplane array
   * @param p0 point that is inside
   * @param vs hyperplane array, normal should go outside
   */
	LinearConstraint(const Vecf<Dim> p0, const vec_E<Hyperplane<Dim>>& vs) {
		const unsigned int size = vs.size();
		MatDNf A(size, Dim);
		VecDf b(size);

		for (unsigned int i = 0; i < size; i++) {
			const auto n = vs[i].n;
			decimal_t c = vs[i].p.dot(n);
			if (n.dot(p0) - c > 0) {
				n = -n;
				c = -c;
			}
			A.row(i) = n;
			b(i) = c;
		}

		A_ = A;
		b_ = b;
	}

  /// Check if the point is inside polyhedron using linear constraint
  bool inside(const Vecf<Dim> &pt) {
    VecDf d = A_ * pt - b_;
    for (unsigned int i = 0; i < d.rows(); i++) {
      if (d(i) > 0)
        return false;
    }
    return true;
  }

  MatDNf<Dim> A_;
  VecDf b_;
};

///LinearConstraint 2D
typedef LinearConstraint<2> LinearConstraint2D;
///LinearConstraint 3D
typedef LinearConstraint<3> LinearConstraint3D;


