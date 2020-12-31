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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <
  typename h_t,
  typename g_t,
  typename r_t,
  typename j_t,
  typename ud_ops_type = void
  >
class HessianGradientOperatorsRJApiNoWeighting
{
public:
  using sc_t  = typename ::pressio::containers::details::traits<h_t>::scalar_t;

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  mutable r_t r_;
  mutable j_t J_;
  g_t g_;
  h_t H_;
  ::pressio::utils::instance_or_reference_wrapper<ud_ops_type> udOps_;

public:
  HessianGradientOperatorsRJApiNoWeighting() = delete;
  HessianGradientOperatorsRJApiNoWeighting(HessianGradientOperatorsRJApiNoWeighting const &) = default;
  HessianGradientOperatorsRJApiNoWeighting & operator=(HessianGradientOperatorsRJApiNoWeighting const &) = default;
  HessianGradientOperatorsRJApiNoWeighting(HessianGradientOperatorsRJApiNoWeighting && o) = default;
  HessianGradientOperatorsRJApiNoWeighting & operator=(HessianGradientOperatorsRJApiNoWeighting && o) = default;
  ~HessianGradientOperatorsRJApiNoWeighting() = default;

  template <
    typename system_t,
    typename state_t,
    typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>,
    mpl::enable_if_t<
      (pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
       pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value)
      and std::is_void<_ud_ops_type>::value,
      int
      > = 0
    >
  HessianGradientOperatorsRJApiNoWeighting(const system_t & system,
					   const state_t & state)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      g_(state),
      H_(::pressio::ops::product<h_t>(pT, pnT,
				      ::pressio::utils::constants<sc_t>::one(),
				      J_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

  template <
    typename system_t,
    typename state_t,
    typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>,
    mpl::enable_if_t<
      (pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
       pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value)
      and !std::is_void<_ud_ops_type>::value,
      int
      > = 0
    >
  HessianGradientOperatorsRJApiNoWeighting(const system_t & system,
					   const state_t & state,
					   _ud_ops_type && udOps)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      g_(state),
      H_(udOps.template product<h_t>(pT, pnT,
				     utils::constants<sc_t>::one(),
				     *J_.data(), *J_.data())),
      udOps_(std::forward<_ud_ops_type>(udOps))
  {}

public:
  void resetForNewCall()    { /* no op */ }
  h_t & hessianRef()     { return H_; }
  g_t & gradientRef()    { return g_; }
  const h_t & hessianCRef() const  { return H_; }
  const g_t & gradientCRef() const { return g_; }

  sc_t getParameter(std::string key) const {
    throw std::runtime_error("GN HessGrad operators does not have parameters");
    return {};
  }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("GN HessGrad operators do not have parameters");
  }

public:
  template<typename system_t, typename state_t>
  mpl::enable_if_t<
  pressio::solvers::constraints::system_residual_jacobian<system_t>::value
  >
  computeOperators(const system_t & systemObj,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    // compute r_
    systemObj.residual(state, r_);

    // compute norm of r_
    residualNorm = this->_computeNormR();

    // recompute Jacobian is needed
    if (recomputeSystemJacobian){
      systemObj.jacobian(state, J_);
      this->_computeHessian();
    }

    // gradient always computed because residual always changes
    this->_computeGradient();
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & systemObj,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    systemObj.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);

    // compute  norm of r_
    residualNorm = this->_computeNormR();

    // hessian only recomputed if Jacobian has been updated
    if (recomputeSystemJacobian){
      this->_computeHessian();
    }

    // gradient always computed because residual always changes
    this->_computeGradient();
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & systemObj,
	       const state_t & state,
	       sc_t & residualNorm) const
  {
    systemObj.residual(state, r_);
    residualNorm = this->_computeNormR();
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & systemObj,
	       const state_t & state,
	       sc_t & residualNorm) const
  {
    systemObj.residualAndJacobian(state, r_, J_, false);
    residualNorm = this->_computeNormR();
  }

private:
  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value, sc_t>
  _computeNormR() const
  {
    return ::pressio::ops::norm2(r_);
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value, sc_t >
  _computeNormR() const
  {
    return udOps_.get().norm2(*r_.data());
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value >
  _computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, beta, H_);
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value >
  _computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    udOps_.get().product(pT, pnT, alpha, *J_.data(), *J_.data(), beta, H_);
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value >
  _computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T r)
    ::pressio::ops::product(pT, alpha, J_, r_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value >
  _computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T r)
    udOps_.get().product(pT, alpha, *J_.data(), *r_.data(), beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_RJ_API_HPP_
