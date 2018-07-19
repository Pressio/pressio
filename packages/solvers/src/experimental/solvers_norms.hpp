
#ifndef SOLVERS_EXPERIMENTAL_NORMS_HPP
#define SOLVERS_EXPERIMENTAL_NORMS_HPP

#include "vector/core_vector_traits_exp.hpp"

namespace solvers {

struct L2Norm {

  template <
    typename T, 
    typename U
  >
  static double compute_norm_difference(const T& lVec, const U& rVec) {

  	typedef typename core::details::vector_traits<T> t_vector_traits_type;
  	typedef typename core::details::vector_traits<U> u_vector_traits_type;

  	static_assert(t_vector_traits_type::vector_class != core::details::WrappedClass::Undefined, "Error: the first argument is not a valid core vector");
  	static_assert(u_vector_traits_type::vector_class != core::details::WrappedClass::Undefined, "Error: the second argument is not a valid core vector");
  	static_assert(core::details::same_vector_structure<T, U>::value, "Error, the two vectors have incompatible dimensions");

    double value = 0;
  	for (auto i = 0; i < lVec.size(); i++) {
  		double tmp = lVec[i] - rVec[i];
  		value += tmp*tmp;
  	}
  	return std::sqrt(value);
  }	
};

} // end namespace solvers

#endif