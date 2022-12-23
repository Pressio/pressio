/*
//@HEADER
// ************************************************************************
//
// solvers_lm_hessian_gradient_operators_rj_api.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_LM_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_LM_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <
  class HessianType,
  class GradientType,
  class ResidualType,
  class JacobianType,
  template<class ...> class hgRJApi_t,
  class ...Args
  >
class LMHessianGradientOperatorsRJApi
{
public:
  using residual_norm_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  // HGOpRJApi_ computes H = J^T J, and g = J^T r
  hgRJApi_t<HessianType, GradientType, ResidualType, JacobianType, Args...>  HGOpRJApi_;

  // lmH = H + lm*diag(H)
  HessianType lmH_;

  // damping factor for LM: we use the same scalar_type of the hessian since
  // the damp factor is used to modify the diagonal of the hessian
  using hessian_scalar_type = typename ::pressio::Traits<HessianType>::scalar_type;
  hessian_scalar_type dampParam_ = pressio::utils::Constants<hessian_scalar_type>::one();

public:
  LMHessianGradientOperatorsRJApi() = delete;
  LMHessianGradientOperatorsRJApi(LMHessianGradientOperatorsRJApi const &) = default;
  LMHessianGradientOperatorsRJApi & operator=(LMHessianGradientOperatorsRJApi const &) = default;
  LMHessianGradientOperatorsRJApi(LMHessianGradientOperatorsRJApi && o) = default;
  LMHessianGradientOperatorsRJApi & operator=(LMHessianGradientOperatorsRJApi && o) = default;
  ~LMHessianGradientOperatorsRJApi() = default;

  template <typename SystemType>
// #ifdef PRESSIO_ENABLE_CXX20
//   requires SystemWithResidualAndJacobian<SystemType> || SystemWithFusedResidualAndJacobian<SystemType> 
// #else
//   template <
//     typename SystemType,
//     mpl::enable_if_t<
//       (::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value or
//        ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value),
//       int
//       > = 0
//     >
// #endif
  LMHessianGradientOperatorsRJApi(const SystemType & system)
    : HGOpRJApi_(system),
      lmH_(::pressio::ops::clone(HGOpRJApi_.hessianCRef()))
  {
    ::pressio::ops::set_zero(lmH_);
  }

  template <
    typename SystemType,
    typename ...ArgsIn,
    mpl::enable_if_t<
      /*(::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value or
       ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value) and */
      sizeof ...(ArgsIn) >= 1,
      int
      > = 0
    >
  LMHessianGradientOperatorsRJApi(const SystemType & system, ArgsIn && ...args)
    : HGOpRJApi_(system, std::forward<ArgsIn>(args)...),
      lmH_(::pressio::ops::clone(HGOpRJApi_.hessianCRef()))
  {
    ::pressio::ops::set_zero(lmH_);
  }

public:
  HessianType & hessianRef()		   { return lmH_; }
  GradientType & gradientRef()		   { return HGOpRJApi_.gradientRef(); }
  const HessianType & hessianCRef() const  { return lmH_; }
  const GradientType & gradientCRef() const { return HGOpRJApi_.gradientCRef(); }

  const HessianType & hessianCRefBeforeLMDiagonalScaling() const {
    return HGOpRJApi_.hessianCRef();
  }

  void setLMDampParam(hessian_scalar_type parIn){ dampParam_ = parIn; }
  auto lmDampParam() const{ return dampParam_; }

public:
  void resetForNewCall(){
    dampParam_ = pressio::utils::Constants<hessian_scalar_type>::one();
  }

  template<typename SystemType, typename StateType>
  void computeOperators(const SystemType & sys,
			const StateType & state,
			residual_norm_type & residualNorm,
			bool recomputeSystemJacobian=true)
  {
    HGOpRJApi_.computeOperators(sys, state,
				residualNorm,
				recomputeSystemJacobian);

    if(recomputeSystemJacobian){
      // compute lmH = H + mu*diag(H)
      const auto & H = HGOpRJApi_.hessianCRef();
      ::pressio::ops::deep_copy(lmH_, H);

      const auto diagH   = ::pressio::diag(H);
      auto diaglmH = ::pressio::diag(lmH_);
      constexpr auto one  = pressio::utils::Constants<hessian_scalar_type>::one();
      ::pressio::ops::update(diaglmH, one, diagH, dampParam_);
    }
  }

  template< typename SystemType, typename StateType>
  void residualNorm(const SystemType & system,
		    const StateType & state,
		    residual_norm_type & residualNorm) const
  {
    HGOpRJApi_.residualNorm(system, state, residualNorm);
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_LM_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
