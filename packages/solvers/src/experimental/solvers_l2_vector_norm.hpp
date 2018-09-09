
#ifndef SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP
#define SOLVERS_EXPERIMENTAL_VECTOR_NORM_HPP

#include "vector/core_vector_traits_exp.hpp"
#include "solvers_meta_static_checks.hpp"


namespace solvers {
	
struct L2Norm {

  template <
    typename T, 
    typename U
  >
  static double compute_norm_difference(const T& lVec, const U& rVec) {

  	typedef typename core::details::vector_traits<T> t_vector_traits_type;
  	typedef typename core::details::vector_traits<U> u_vector_traits_type;

  	static_assert(t_vector_traits_type::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined, "Error: the first argument is not a valid core vector");
  	static_assert(u_vector_traits_type::wrapped_package_identifier != core::details::WrappedPackageIdentifier::Undefined, "Error: the second argument is not a valid core vector");
  	static_assert(solvers::meta::same_vector_structure<T, U>::value, "Error, the two vectors have incompatible dimensions");

    double value = 0;

    auto dVec = lVec - rVec;
    dVec.norm2(value);
  	
  	return value;
  }	
};


} // end namespace solvers

#endif
