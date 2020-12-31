/*
//@HEADER
// ************************************************************************
//
// solvers_gn_hessian_gradient_operators_weighted_rj_api.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WEIGHTED_RJ_API_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WEIGHTED_RJ_API_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <class r_t, class j_t, class T, class operand_t, class result_t>
mpl::enable_if_t<
  ::pressio::solvers::constraints::least_squares_weighting_operator_accepting_wrappers<T,r_t,j_t>::value
  >
_applyWeightingHelper(const T & functorM,
		      const operand_t & operand,
		      result_t & result,
		      bool is_irwls,
		      int callCount)
{
  if(is_irwls && callCount > 1)
  {
    functorM(operand, result);
  }else{
    functorM(operand, result);
  }
};

template <class r_t, class j_t, class T, class operand_t, class result_t>
mpl::enable_if_t<
  ::pressio::solvers::constraints::least_squares_weighting_operator_accepting_native<T,r_t,j_t>::value
  >
_applyWeightingHelper(const T & functorM,
		      const operand_t & operand,
		      result_t & result,
		      bool is_irwls,
		      int callCount)
{
  if(is_irwls && callCount > 1)
  {
    functorM(*operand.data(), *result.data());
  }else{
    functorM(*operand.data(), *result.data());
  }
};

template <class r_t, class j_t, class T, class operand_t, class result_t>
mpl::enable_if_t<
  ::pressio::solvers::constraints::least_squares_weighting_operator_accepting_wrappers<T,r_t,j_t>::value
  >
_applyWeightingHelper(const T & functorM,
		      const operand_t & operand,
		      result_t & result)
{
  functorM(operand, result);
};

template <class r_t, class j_t, class T, class operand_t, class result_t>
mpl::enable_if_t<
  ::pressio::solvers::constraints::least_squares_weighting_operator_accepting_native<T,r_t,j_t>::value
  >
_applyWeightingHelper(const T & functorM,
		      const operand_t & operand,
		      result_t & result)
{
  functorM(*operand.data(), *result.data());
};




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

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  int callCount_ = 0;
  mutable r_t r_;
  mutable r_t Mr_;
  j_t J_;
  j_t MJ_;
  g_t g_;
  h_t H_;
  ::pressio::utils::instance_or_reference_wrapper<ud_ops_type> udOps_;
  ::pressio::utils::instance_or_reference_wrapper<weighting_functor_t> functorM_;

  static constexpr auto is_irwls =
    std::is_same<
    weighting_functor_t,
    ::pressio::solvers::nonlinear::impl::IrwWeightingOperator<r_t, j_t>
    >::value;

public:
  WeightedHessianGradientOperatorsRJApi() = delete;
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi const &) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi const &) = default;
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi && o) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi && o) = default;
  ~WeightedHessianGradientOperatorsRJApi() = default;

  template <
    typename system_t,
    typename state_t,
    typename _weigh_t,
    typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>,
    mpl::enable_if_t<
      (pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
       pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value)
      and std::is_void<_ud_ops_type>::value,
      int
      > = 0
    >
  WeightedHessianGradientOperatorsRJApi(const system_t & system,
					const state_t & state,
					_weigh_t && functorM)
  : r_(system.createResidual()),
    Mr_(r_),
    J_(system.createJacobian()),
    MJ_(J_),
    g_(state),
    H_(::pressio::ops::product<h_t>(pT, pnT,
				    ::pressio::utils::constants<sc_t>::one(),
				    J_)),
    functorM_(std::forward<_weigh_t>(functorM))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(Mr_);
    ::pressio::ops::set_zero(J_);
    ::pressio::ops::set_zero(MJ_);
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

  template <
    typename system_t,
    typename state_t,
    typename _weigh_t,
    typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>,
    mpl::enable_if_t<
      (pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
       pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value)
      and !std::is_void<_ud_ops_type>::value,
      int
      > = 0
    >
  WeightedHessianGradientOperatorsRJApi(const system_t & system,
					const state_t & state,
					_ud_ops_type && udOps,
					_weigh_t && functorM)
  : r_(system.createResidual()),
    Mr_(r_),
    J_(system.createJacobian()),
    MJ_(J_),
    g_(state),
    H_(udOps.template product<h_t>(pT, pnT,
				   utils::constants<sc_t>::one(),
				   *J_.data(), *J_.data())),
    udOps_(std::forward<_ud_ops_type>(udOps)),
    functorM_(std::forward<_weigh_t>(functorM))
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <
    typename system_t,
    typename state_t,
    typename _weigh_t = weighting_functor_t,
    typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>,
    mpl::enable_if_t<
      (pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
       pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value)
      and std::is_void<_ud_ops_type>::value
      and std::is_same<_weigh_t,
		       ::pressio::solvers::nonlinear::impl::IrwWeightingOperator<r_t, j_t>
		       >::value
      , int > = 0
    >
  WeightedHessianGradientOperatorsRJApi(const system_t & system,
					const state_t & state,
					sc_t pValue)
  : r_(system.createResidual()),
    Mr_(r_),
    J_(system.createJacobian()),
    MJ_(J_),
    g_(state),
    H_(::pressio::ops::product<h_t>(pT, pnT,
				    ::pressio::utils::constants<sc_t>::one(),
				    J_)),
    functorM_(std::move(_weigh_t(system)))
  {
    this->set_p(pValue);
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(Mr_);
    ::pressio::ops::set_zero(J_);
    ::pressio::ops::set_zero(MJ_);
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }
#endif

public:
  void resetForNewCall() {
    callCount_ = 0;
  }

  h_t & hessianRef()		   { return H_; }
  g_t & gradientRef()		   { return g_; }
  const h_t & hessianCRef() const  { return H_; }
  const g_t & gradientCRef() const { return g_; }

  template<bool _is_irwls = is_irwls>
  mpl::enable_if_t<_is_irwls> set_p(sc_t pIn)
  {
    functorM_.get().set_p(pIn);
  }

  sc_t getParameter(std::string key) const
  {
    throw std::runtime_error("GN HessGrad operators does not have parameters");
    return {};
  }

  template <typename T>
  void setParameter(std::string key, T value)
  {
    throw std::runtime_error("GN HessGrad operators do not have parameters");
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & system,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    callCount_++;

    // compute r from system object
    system.residual(state, r_);
    // apply M
    _applyWeightingHelper<r_t,j_t>(functorM_.get(), r_, Mr_, is_irwls, callCount_);

    residualNorm = this->_computeNorm();

    if (recomputeSystemJacobian){
      system.jacobian(state, J_);
      _applyWeightingHelper<r_t,j_t>(functorM_.get(), J_, MJ_, is_irwls, callCount_);
      _computeHessian();
    }

    _computeGradient();
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & system,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    callCount_++;

    system.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    _applyWeightingHelper<r_t,j_t>(functorM_.get(), r_, Mr_, is_irwls, callCount_);
    residualNorm = this->_computeNorm();

    if (recomputeSystemJacobian){
      _applyWeightingHelper<r_t,j_t>(functorM_.get(), J_, MJ_, is_irwls, callCount_);
      _computeHessian();
    }

    _computeGradient();
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system,
	       const state_t & state,
	       sc_t & residualNorm) const
  {
    system.residual(state, r_);
    _applyWeightingHelper<r_t,j_t>(functorM_.get(), r_, Mr_);
    residualNorm = this->_computeNorm();
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system,
	       const state_t & state,
	       sc_t & residualNorm) const
  {
    // here we query system to recompute r_ only (that is why we pass false)
    system.residualAndJacobian(state, r_, J_, false);
    _applyWeightingHelper<r_t,j_t>(functorM_.get(), r_, Mr_);
    residualNorm = this->_computeNorm();
  }

private:
  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value, sc_t >
  _computeNorm() const
  {
    return std::sqrt(::pressio::ops::dot(r_, Mr_));
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value, sc_t >
  _computeNorm() const
  {
    return std::sqrt(udOps_.get().dot(r_, Mr_));
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value >
  _computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, MJ_, beta, H_);
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value >
  _computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    udOps_.get().product(pT, pnT, alpha, *J_.data(), *MJ_.data(), beta, H_);
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< std::is_void<_ud_ops_type>::value >
  _computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T M r)
    ::pressio::ops::product(pT, alpha, J_, Mr_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename _ud_ops_type = mpl::remove_cvref_t<ud_ops_type>>
  mpl::enable_if_t< mpl::not_void<_ud_ops_type>::value >
  _computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::constants<sc_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<sc_t>::one();
    // compute gradient (g_ = J^T M r)
    udOps_.get().product(pT, alpha, *J_.data(), *Mr_.data(), beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WEIGHTED_RJ_API_HPP_
