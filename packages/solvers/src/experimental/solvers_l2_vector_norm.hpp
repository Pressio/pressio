
#ifndef SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP
#define SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP

#include "solvers_meta_static_checks.hpp"
#include "../../../algebra/src/containers_ops/norms/algebra_norm2_vector.hpp"


namespace rompp{
namespace solvers{

struct L2Norm {

  template <
    typename T,
    typename U
  >
  static double compute_norm_difference(const T& lVec, const U& rVec) {

    static_assert(solvers::meta::are_vector_compatible<T, U>::value, "Error: the two vectors are not compatible");

    T dVec(lVec - rVec);
    double value =::rompp::algebra::ops::norm2(dVec);
    // dVec.norm2(value);
    return value;
  }


	template <
    typename T
  >
  static double compute_norm(const T& vec) {
		double value = ::rompp::algebra::ops::norm2(vec);
    // vec.norm2(value);
		return value;
  }
};


} // end namespace solvers

}//end namespace rompp
#endif
