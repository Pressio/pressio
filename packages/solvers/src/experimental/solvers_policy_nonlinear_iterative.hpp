#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_ITERATIVE_HPP

#include <iostream>
#include <type_traits>

#include "system_traits.hpp"
#include "meta/core_meta_static_checks.hpp"


namespace solvers {

struct SolversNonLinearIterativeNewtonRaphsonPolicy {
  

  template <
    typename SystemT,
    typename VectorT,
    typename std::enable_if<
      !core::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value,
      int
    >::type* = nullptr
  >
  static auto solve(const SystemT& system, const VectorT& b) {

    std::cerr << "Error: the type of the RHS vector is not compatible with the provided nonlinear system" << std::endl;
    assert(0);

  	return b;
  }


  template <
    typename SystemT,
    typename VectorT,
    typename std::enable_if<
      core::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value,
      int
    >::type* = nullptr
  >
  static auto solve(const SystemT& sytem, const VectorT& b) {
  	return 0;
  }
};

} // end namespace solvers

#endif