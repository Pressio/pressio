/*
//@HEADER
// ************************************************************************
//
// solvers_gn_hessian_gradient_operators_with_weighting.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WITH_WEIGHTING_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WITH_WEIGHTING_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{


template <
  typename h_t, 
  typename g_t,
  typename r_t, 
  typename j_t,
  typename ud_ops_type = void,
  typename weighting_functor_t = void
  >
class WeightedHessianGradientOperatorsRJApi
{
public:
  using sc_t  = typename ::pressio::containers::details::traits<h_t>::scalar_t;
  using ud_ops_t = ud_ops_type;

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  mutable r_t r_;
  mutable r_t Mr_;

  j_t J_;
  j_t MJ_;

  g_t g_;
  h_t H_;

  const ud_ops_t * udOps_   = nullptr;
  const weighting_functor_t * functorM_= nullptr;

public:
  WeightedHessianGradientOperatorsRJApi() = delete;

  template <
   typename system_t, 
   typename state_t, 
   typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
       pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value)
      and std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  WeightedHessianGradientOperatorsRJApi(const system_t & system,
        const state_t & state,
        const weighting_functor_t & functorM)
    : r_(system.createResidual()),
      Mr_(r_),
      J_(system.createJacobian()),
      MJ_(J_),
      g_(state),
      H_(::pressio::ops::product<h_t>(pT, pnT,
              ::pressio::utils::constants<sc_t>::one(),
              J_)),
      functorM_(&functorM)
  {}

  template <
   typename system_t, 
   typename state_t, 
   typename _ud_ops_t = ud_ops_t,
    mpl::enable_if_t<
      (pressio::solvers::concepts::system_residual_jacobian<system_t>::value or
       pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value)
      and !std::is_void<_ud_ops_t>::value,
      int
     > = 0
  >
  WeightedHessianGradientOperatorsRJApi(const system_t & system,
        const state_t & state,
        const _ud_ops_t & udOps,
        const weighting_functor_t & functorM)
    : r_(system.createResidual()),
      Mr_(r_),
      J_(system.createJacobian()),
      MJ_(J_),
      g_(state),
      H_(udOps.template product<h_t>(pT, pnT,
             utils::constants<sc_t>::one(),
             *J_.data(), *J_.data())),
      udOps_(&udOps),
      functorM_(&functorM)
  {}

  // copy constr and assign
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi const &) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi const &) = default;

  // move constr and assign
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi && o) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi && o) = default;

  // destr
  ~WeightedHessianGradientOperatorsRJApi() = default;

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

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & system,
       const state_t & state,
       ::pressio::Norm normType,
       sc_t & residualNorm,
       bool recomputeSystemJacobian = true)
  {
    // compute r from system object
    system.residual(state, r_);//, normType, residualNorm);
    // apply M 
    (*functorM_)(r_, Mr_);

    assert(normType == ::pressio::Norm::L2);
    residualNorm = this->computeNorm(normType);

    if (recomputeSystemJacobian){
      system.jacobian(state, J_);
      (*functorM_)(J_, MJ_);
      computeHessian();
    }

    computeGradient();
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & system,
       const state_t & state,
       ::pressio::Norm normType,
       sc_t & residualNorm,
       bool recomputeSystemJacobian = true)
  {
    system.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    //normType, residualNorm, recomputeSystemJacobian);

    (*functorM_)(r_, Mr_);
    assert(normType == ::pressio::Norm::L2);
    residualNorm = this->computeNorm(normType);

    if (recomputeSystemJacobian){
      (*functorM_)(J_, MJ_);
      computeHessian();
    }

    computeGradient();
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system, 
         const state_t & state,
         ::pressio::Norm normType, 
         sc_t & residualNorm) const
  {
    system.residual(state, r_);//, normType, residualNorm);
    (*functorM_)(r_, Mr_);
    residualNorm = this->computeNorm(normType);
    // residualNorm = ::pressio::ops::dot(r_, Mr_);
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::concepts::system_fused_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system, 
         const state_t & state,
         ::pressio::Norm normType, 
         sc_t & residualNorm) const
  {
    // system.residualNorm(state, normType, residualNorm);

    // here we query system to recompute r_ only (that is why we pass false)
    system.residualAndJacobian(state, r_, J_, false); 
                              //normType,residualNorm, false);
    (*functorM_)(r_, Mr_);
    residualNorm = this->computeNorm(normType);
    // residualNorm = ::pressio::ops::dot(r_, Mr_);
  }

private:
  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value, sc_t >
  computeNorm(::pressio::Norm normType) const
  {
    return ::pressio::ops::dot(r_, Mr_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<_ud_ops_t>::value, sc_t >
  computeNorm(::pressio::Norm normType) const
  {
    return udOps_->dot(r_, Mr_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, MJ_, beta, H_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T M r)
    ::pressio::ops::product(pT, alpha, J_, Mr_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    udOps_->product(pT, pnT, alpha, *J_.data(), *MJ_.data(), beta, H_);
  }

  template<typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T M r)
    udOps_->product(pT, alpha, *J_.data(), *Mr_.data(), beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WITH_WEIGHTING_HPP_
