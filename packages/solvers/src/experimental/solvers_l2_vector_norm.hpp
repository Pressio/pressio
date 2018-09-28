
#ifndef SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP
#define SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP

#include "solvers_meta_static_checks.hpp"


namespace rompp{
namespace solvers{

struct L2Norm {

  template <
    typename T,
    typename U
  >
  static double compute_norm_difference(const T& lVec, const U& rVec) {

    static_assert(solvers::meta::are_vector_compatible<T, U>::value, "Error: the two vectors are not compatible");

    double value = 0;
    T dVec(lVec - rVec);

    dVec.norm2(value);
    return value;
  }


	template <
    typename T
  >
  static double compute_norm(const T& vec) {
		double value;
    vec.norm2(value);
		return value;
  }
};


} // end namespace solvers

}//end namespace rompp
#endif
