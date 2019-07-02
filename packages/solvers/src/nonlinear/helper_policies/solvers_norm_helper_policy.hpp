
#ifndef SOLVERS_IMPL_NORM_HELPER_POLICY_HPP
#define SOLVERS_IMPL_NORM_HELPER_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CONTAINERS_OPS"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename convergence_tag>
struct NormSelectorHelper{
  using norm_t = L2Norm;
};

template <typename norm_type>
struct NormSelectorHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_type>
  >{
  using norm_t = norm_type;
};
//---------------------------------------------------------


template <typename norm_t>
struct ComputeNormHelper;

template <>
struct ComputeNormHelper<::pressio::solvers::L2Norm>{
  template <typename vec_t, typename scalar_t>
  static void evaluate(const vec_t & vecIn, scalar_t & result) {
    result = ::pressio::containers::ops::norm2(vecIn);
  }
};

template <>
struct ComputeNormHelper<::pressio::solvers::L1Norm>{
  template <typename vec_t, typename scalar_t>
  static void evaluate(const vec_t & vecIn, scalar_t & result) {
    result = ::pressio::containers::ops::norm1(vecIn);
  }
};
//---------------------------------------------------------


}}}} //end namespace pressio::solvers::iterative::impl
#endif
