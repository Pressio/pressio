/*
//@HEADER
// ************************************************************************
//
// solvers_gn_hessian_gradient_operators_rj_api.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_FUSED_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_FUSED_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class HessianType, class GradientType, class ResidualType, class JacobianType>
class FusedHessianGradientOperatorsRJApiNoWeighting
{
public:
  using residual_norm_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  mutable ResidualType r_;
  mutable JacobianType J_;
  GradientType g_;
  HessianType H_;

public:
  FusedHessianGradientOperatorsRJApiNoWeighting() = delete;
  FusedHessianGradientOperatorsRJApiNoWeighting(FusedHessianGradientOperatorsRJApiNoWeighting const &) = default;
  FusedHessianGradientOperatorsRJApiNoWeighting & operator=(FusedHessianGradientOperatorsRJApiNoWeighting const &) = default;
  FusedHessianGradientOperatorsRJApiNoWeighting(FusedHessianGradientOperatorsRJApiNoWeighting && o) = default;
  FusedHessianGradientOperatorsRJApiNoWeighting & operator=(FusedHessianGradientOperatorsRJApiNoWeighting && o) = default;
  ~FusedHessianGradientOperatorsRJApiNoWeighting() = default;

  template <typename SystemType>
  FusedHessianGradientOperatorsRJApiNoWeighting(const SystemType & system)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      g_(system.createState()),
      H_(::pressio::ops::product<HessianType>
	 (pT, pnT, ::pressio::utils::Constants<
	  typename ::pressio::Traits<JacobianType>::scalar_type>::one(),
	  J_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

public:
  void resetForNewCall()    { /* no op */ }
  HessianType & hessianRef()     { return H_; }
  GradientType & gradientRef()    { return g_; }
  const HessianType & hessianCRef() const  { return H_; }
  const GradientType & gradientCRef() const { return g_; }

  residual_norm_type getParameter(std::string key) const {
    throw std::runtime_error("GN HessGrad operators does not have parameters");
    return {};
  }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("GN HessGrad operators do not have parameters");
  }

public:

  template<typename SystemType, typename StateType>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithFusedResidualAndJacobian<SystemType>
  void 
#else  
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value
  >  
#endif
  computeOperators(const SystemType & systemObj,
		   const StateType & state,
		   residual_norm_type & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    systemObj.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);

    // compute  norm of r_
    residualNorm = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }

    // hessian only recomputed if Jacobian has been updated
    if (recomputeSystemJacobian){
      this->_computeHessian();
    }

    // gradient always computed because residual always changes
    this->_computeGradient();
  }

  template< typename SystemType, typename StateType>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithResidualAndJacobian<SystemType>
  void 
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value
  >  
#endif
  residualNorm(const SystemType & systemObj,
	       const StateType & state,
	       residual_norm_type & residualNorm) const
  {
    systemObj.residual(state, r_);
    residualNorm = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

  template< typename SystemType, typename StateType>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithFusedResidualAndJacobian<SystemType>
  void
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value
  >  
#endif
  residualNorm(const SystemType & systemObj,
	       const StateType & state,
	       residual_norm_type & residualNorm) const
  {
    systemObj.residualAndJacobian(state, r_, J_, false);
    residualNorm = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

private:
  void _computeHessian()
  {
    using H_scalar_type = typename ::pressio::Traits<HessianType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<H_scalar_type>::zero();
    using J_scalar_type = typename ::pressio::Traits<JacobianType>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, beta, H_);
  }

  void _computeGradient()
  {
    using g_scalar_type = typename ::pressio::Traits<GradientType>::scalar_type;
    constexpr auto beta  = ::pressio::utils::Constants<g_scalar_type>::zero();
    using J_scalar_type = typename ::pressio::Traits<JacobianType>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<J_scalar_type>::one();
    // compute gradient (g_ = J^T r)
    ::pressio::ops::product(pT, alpha, J_, r_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<g_scalar_type>::negOne());
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
