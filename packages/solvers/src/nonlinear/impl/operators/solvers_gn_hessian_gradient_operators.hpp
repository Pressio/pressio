/*
//@HEADER
// ************************************************************************
//
// solvers_hessian_gradient_operators.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename h_t, typename g_t>
class HessianGradientOperatorsHGApi
{
  using sc_t = typename ::pressio::containers::details::traits<h_t>::scalar_t;
  g_t g_;
  h_t H_;

public:
  HessianGradientOperatorsHGApi() = delete;

  template <
   typename system_t, typename state_t,
    mpl::enable_if_t<
      pressio::solvers::concepts::system_hessian_gradient<system_t>::value or
      pressio::solvers::concepts::system_fused_hessian_gradient<system_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsHGApi(const system_t & system, const state_t & state)
    : g_( system.createGradient() ),
      H_( system.createHessian() )
  {}

public:
  h_t & getHessian(){ return H_; }
  g_t & getGradient(){ return g_; }
  const h_t & getHessian() const { return H_; }
  const g_t & getGradient() const{ return g_; }

  void resetForNewCall(){
    //no op
  }

  template< typename system_t, typename state_t>
  void residualNorm(const system_t & system, const state_t & state,
		    ::pressio::Norm normType, sc_t & residualNorm) const
  {
    system.residualNorm(state, normType, residualNorm);
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_hessian_gradient<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   ::pressio::Norm normType,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    if (recomputeSystemJacobian){
      sys.hessian(state, H_);
    }

    sys.gradient(state, g_, normType,
		 residualNorm,
		 recomputeSystemJacobian);

    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_fused_hessian_gradient<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   ::pressio::Norm normType,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    sys.hessianAndGradient(state, H_, g_, normType,
			   residualNorm, recomputeSystemJacobian);

    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};



template <
  typename h_t, typename g_t,
  typename r_t, typename j_t,
  typename ud_ops_t = void
  >
class HessianGradientOperatorsRJApi
{
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();
  using sc_t = typename ::pressio::containers::details::traits<h_t>::scalar_t;

  mutable r_t r_;
  j_t J_;
  g_t g_;
  h_t H_;
  const ud_ops_t * udOps_ = nullptr;

public:
  HessianGradientOperatorsRJApi() = delete;

  template <
    typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
       pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value)
      and std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsRJApi(const system_t & system, const state_t & state)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      g_(state),
      H_(::pressio::ops::product<h_t>(pT, pnT,
				      ::pressio::utils::constants<sc_t>::one(),
				      J_))
  {}

  template <
   typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
       pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value)
      and !std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsRJApi(const system_t & system,
				const state_t & state,
				const _ud_ops_t & udOps)
    : r_(system.createResidual()),
      J_(system.createJacobian()),
      g_(state),
      H_(udOps.template product<h_t>(pT, pnT,
				     utils::constants<sc_t>::one(),
				     *J_.data(), *J_.data())),
      udOps_(&udOps)
  {}

public:
  h_t & getHessian(){ return H_; }
  g_t & getGradient(){ return g_; }
  const h_t & getHessian() const { return H_; }
  const g_t & getGradient() const{ return g_; }

  void resetForNewCall(){
    //no op
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system, const state_t & state,
	       ::pressio::Norm normType, sc_t & residualNorm) const
  {
    system.residual(state, r_, normType, residualNorm);
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system, const state_t & state,
	       ::pressio::Norm normType, sc_t & residualNorm)
  {
    system.residualNorm(state, normType, residualNorm);
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   ::pressio::Norm normType,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    sys.residual(state, r_, normType, residualNorm);

    if (recomputeSystemJacobian){
      sys.jacobian(state, J_);
      computeHessian();
    }

    computeGradient();
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   ::pressio::Norm normType,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    sys.residualAndJacobian(state, r_, J_, normType,
			    residualNorm,
			    recomputeSystemJacobian);
    if (recomputeSystemJacobian){
      computeHessian();
    }

    computeGradient();
  }

private:
  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  computeHessian(){
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, beta, H_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  computeGradient(){
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T r)
    ::pressio::ops::product(pT, alpha, J_, r_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  computeHessian(){
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    udOps_->product(pT, pnT, alpha, *J_.data(), *J_.data(), beta, H_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  computeGradient(){
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T r)
    udOps_->product(pT, alpha, *J_.data(), *r_.data(), beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_HESSIAN_GRADIENT_OPERATORS_HPP_
