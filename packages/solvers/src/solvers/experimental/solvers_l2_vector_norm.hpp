
#ifndef SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP
#define SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP

#include "solvers_meta_static_checks.hpp"


namespace solvers {
	
struct L2Norm {

  template <
    typename T, 
    typename U
  >
  static double compute_norm_difference(const T& lVec, const U& rVec) {

    static_assert(solvers::meta::are_vector_compatible<T, U>::value, "Error: the two vectors are not compatible");

    double value = 0;
    auto dVec = lVec - rVec;

    dVec.norm2(value);  	
    return value;
  }	
};


} // end namespace solvers

#endif
