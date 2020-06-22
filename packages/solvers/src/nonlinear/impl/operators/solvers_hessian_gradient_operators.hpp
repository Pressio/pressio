/*
//@HEADER
// ************************************************************************
//
// solvers_hessian_gradient_operators_impl.hpp
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

#ifndef PRESSIO_HESSIAN_GRADIENT_OPERATORS_HPP_
#define PRESSIO_HESSIAN_GRADIENT_OPERATORS_HPP_

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
      pressio::solvers::meta::system_meets_hessian_gradient_api<system_t>::value or
      pressio::solvers::meta::system_meets_fused_hessian_gradient_api<system_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsHGApi(const system_t & system, const state_t & state)
    : g_(system.createGradientObject(state)),
      H_(system.createHessianObject(state)){}

public:
  h_t & getHessian(){ return H_; }
  g_t & getGradient(){ return g_; }
  const h_t & getHessian() const { return H_; }
  const g_t & getGradient() const{ return g_; }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<pressio::solvers::meta::system_meets_hessian_gradient_api<system_t>::value>
  computeOperators(const system_t & sys, const state_t & state,
		   ::pressio::solvers::Norm normType, sc_t & residualNorm)
  {
    sys.hessian(state, H_);
    sys.gradient(state, g_, normType, residualNorm);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<pressio::solvers::meta::system_meets_fused_hessian_gradient_api<system_t>::value>
  computeOperators(const system_t & sys, const state_t & state,
		   ::pressio::solvers::Norm normType, sc_t & residualNorm)
  {
    sys.hessianAndGradient(state, H_, g_, normType, residualNorm);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};



template <typename h_t, typename g_t, typename r_t, typename j_t, typename ud_ops_t = void>
class HessianGradientOperatorsRJApi
{
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();
  using sc_t = typename ::pressio::containers::details::traits<h_t>::scalar_t;

  r_t r_;
  j_t J_;
  g_t g_;
  h_t H_;
  const ud_ops_t * udOps_ = nullptr;

public:
  HessianGradientOperatorsRJApi() = delete;

  template <
    typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
       pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value) and
      std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsRJApi(const system_t & system, const state_t & state)
    : r_(system.createResidualObject(state)),
      J_(system.createJacobianObject(state)),
      g_(state),
      H_(::pressio::ops::product<h_t>(pT, pnT, ::pressio::utils::constants<sc_t>::one(), J_)){}

  template <
   typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
       pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value) and
      !std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsRJApi(const system_t & system, const state_t & state, const ud_ops_t * udOps)
    : r_(system.createResidualObject(state)),
      J_(system.createJacobianObject(state)),
      g_(state),
      H_(udOps.template product<h_t>(pT, pnT, ::pressio::utils::constants<sc_t>::one(), *J_.data(), *J_.data())){}

public:
  h_t & getHessian(){ return H_; }
  g_t & getGradient(){ return g_; }
  const h_t & getHessian() const { return H_; }
  const g_t & getGradient() const{ return g_; }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value>
  computeOperators(const system_t & sys, const state_t & state,
		   ::pressio::solvers::Norm normType, sc_t & residualNorm)
  {
    sys.residual(state, r_, normType, residualNorm);
    sys.jacobian(state, J_);
    computeHessian();
    computeGradient();
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value>
  computeOperators(const system_t & sys, const state_t & state,
		   ::pressio::solvers::Norm normType, sc_t & residualNorm)
  {
    sys.residualAndJacobian(state, r_, J_, normType, residualNorm);
    computeHessian();
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


template <typename h_t, typename g_t, typename r_t, typename j_t, typename ud_ops_t = void>
class LMHessianGradientOperatorsRJApi
{
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();
  using sc_t = typename ::pressio::containers::details::traits<h_t>::scalar_t;

  // HGOpRJApi_ contains H = J^T J, and g = J^T r
  HessianGradientOperatorsRJApi<h_t, g_t, r_t, j_t, ud_ops_t> HGOpRJApi_;

  // lmH contains H + lm*diag(H)
  h_t lmH_;

  // damping factor for LM
  sc_t dampParam_ = pressio::utils::constants<sc_t>::one();

public:
  LMHessianGradientOperatorsRJApi() = delete;

  template <
   typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
       pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value) and
      std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  LMHessianGradientOperatorsRJApi(const system_t & system, const state_t & state)
    : HGOpRJApi_(system, state),
      lmH_(HGOpRJApi_.getHessian()){}

  template <
   typename system_t, typename state_t, typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::meta::system_meets_residual_jacobian_api<system_t>::value or
       pressio::solvers::meta::system_meets_fused_residual_jacobian_api<system_t>::value) and
      !std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  LMHessianGradientOperatorsRJApi(const system_t & system, const state_t & state, const ud_ops_t * udOps)
    : HGOpRJApi_(system, state, udOps),
      lmH_(HGOpRJApi_.getHessian()){}

public:
  h_t & getHessian(){ return lmH_; }
  g_t & getGradient(){ return HGOpRJApi_.getGradient(); }
  const h_t & getHessian() const { return lmH_; }
  const g_t & getGradient() const{ return HGOpRJApi_.getGradient(); }

  const h_t & getHessianBeforeLMDiagonalScaling() const { return HGOpRJApi_.getHessian(); }
  void setLMDampParam(sc_t parIn){ dampParam_ = parIn; }
  sc_t getLMDampParam() const{ return dampParam_; }

public:
  template<typename system_t, typename state_t>
  void computeOperators(const system_t & sys, const state_t & state,
			::pressio::solvers::Norm normType, sc_t & residualNorm)
  {
    HGOpRJApi_.computeOperators(sys, state, normType, residualNorm);

    // compute lmH = H + mu*diag(H)
    const auto & H = HGOpRJApi_.getHessian();
    ::pressio::ops::deep_copy(lmH_, H);
    // this needs to be replaces with a daxpy when we have the diagonal view
    lmH_.data()->diagonal() = lmH_.data()->diagonal() + dampParam_*lmH_.data()->diagonal();
  }
};

}}}}
#endif
