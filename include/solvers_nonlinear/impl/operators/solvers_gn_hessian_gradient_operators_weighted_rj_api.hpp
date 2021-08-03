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

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class ResidualType, class JacobianType, class T, class OperandType, class ResultType>
void _applyWeightingHelper(const T & functorM,
		      const OperandType & operand,
		      ResultType & result,
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

template <class ResidualType, class JacobianType, class T, class OperandType, class ResultType>
void _applyWeightingHelper(const T & functorM, const OperandType & operand, ResultType & result)
{
  functorM(operand, result);
};

template <
  typename HessianType,
  typename GradientType,
  typename ResidualType,
  typename JacobianType,
  typename scalarType,
  typename weighting_functor_t = void
  >
class WeightedHessianGradientOperatorsRJApi
{
public:
  using scalar_type = scalarType;

private:
  static constexpr auto pT  = ::pressio::transpose();
  static constexpr auto pnT = ::pressio::nontranspose();

  int callCount_ = 0;
  mutable ResidualType r_;
  mutable ResidualType Mr_;
  JacobianType J_;
  JacobianType MJ_;
  GradientType g_;
  HessianType H_;
  ::pressio::utils::InstanceOrReferenceWrapper<weighting_functor_t> functorM_;

  static constexpr auto is_irwls =
    std::is_same<
    weighting_functor_t,
    ::pressio::nonlinearsolvers::impl::IrwWeightingOperator<ResidualType, JacobianType, scalar_type>
    >::value;

public:
  WeightedHessianGradientOperatorsRJApi() = delete;
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi const &) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi const &) = default;
  WeightedHessianGradientOperatorsRJApi(WeightedHessianGradientOperatorsRJApi && o) = default;
  WeightedHessianGradientOperatorsRJApi & operator=(WeightedHessianGradientOperatorsRJApi && o) = default;
  ~WeightedHessianGradientOperatorsRJApi() = default;

  template <
    typename SystemType,
    typename StateType,
    typename _weigh_t,
    mpl::enable_if_t<
      (::pressio::nonlinearsolvers::constraints::system_residual_jacobian<SystemType>::value or
       ::pressio::nonlinearsolvers::constraints::system_fused_residual_jacobian<SystemType>::value),
      int
      > = 0
    >
  WeightedHessianGradientOperatorsRJApi(const SystemType & system,
					const StateType & state,
					_weigh_t && functorM)
  : r_(system.createResidual()),
    Mr_(::pressio::ops::clone(r_)),
    J_(system.createJacobian()),
    MJ_(::pressio::ops::clone(J_)),
    g_(::pressio::ops::clone(state)),
    H_(::pressio::ops::product<HessianType>(pT, pnT, ::pressio::utils::Constants<scalar_type>::one(), J_)),
    functorM_(std::forward<_weigh_t>(functorM))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(Mr_);
    ::pressio::ops::set_zero(J_);
    ::pressio::ops::set_zero(MJ_);
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   template <
//     typename SystemType,
//     typename StateType,
//     typename _weigh_t = weighting_functor_t,
//     mpl::enable_if_t<
//       (::pressio::nonlinearsolvers::constraints::system_residual_jacobian<SystemType>::value or
//        ::pressio::nonlinearsolvers::constraints::system_fused_residual_jacobian<SystemType>::value)
//       and std::is_same<_weigh_t,
// 		       ::pressio::nonlinearsolvers::impl::IrwWeightingOperator<ResidualType, JacobianType>
// 		       >::value
//       , int > = 0
//     >
//   WeightedHessianGradientOperatorsRJApi(const SystemType & system,
// 					const StateType & state,
// 					scalar_type pValue)
//   : r_(system.createResidual()),
//     Mr_(::pressio::ops::clone(r_)),
//     J_(system.createJacobian()),
//     MJ_(::pressio::ops::clone(J_)),
//     g_(::pressio::ops::clone(state)),
//     H_(::pressio::ops::product<HessianType>(pT, pnT,::pressio::utils::Constants<scalar_type>::one(), J_)),
//     functorM_(std::move(_weigh_t(system)))
//   {
//     this->set_p(pValue);
//     ::pressio::ops::set_zero(r_);
//     ::pressio::ops::set_zero(Mr_);
//     ::pressio::ops::set_zero(J_);
//     ::pressio::ops::set_zero(MJ_);
//     ::pressio::ops::set_zero(g_);
//     ::pressio::ops::set_zero(H_);
//   }
// #endif

public:
  void resetForNewCall() {
    callCount_ = 0;
  }

  HessianType & hessianRef()		   { return H_; }
  GradientType & gradientRef()		   { return g_; }
  const HessianType & hessianCRef() const  { return H_; }
  const GradientType & gradientCRef() const { return g_; }

  template<bool _is_irwls = is_irwls>
  mpl::enable_if_t<_is_irwls> set_p(scalar_type pIn)
  {
    functorM_.get().set_p(pIn);
  }

  scalar_type getParameter(std::string key) const
  {
    throw std::runtime_error("GN HessGrad operators does not have parameters");
    return {};
  }

  template <typename T>
  void setParameter(std::string key, T value)
  {
    throw std::runtime_error("GN HessGrad operators do not have parameters");
  }

  template<typename SystemType, typename StateType>
  mpl::enable_if_t<
    ::pressio::nonlinearsolvers::constraints::system_residual_jacobian<SystemType>::value
    >
  computeOperators(const SystemType & system,
		   const StateType & state,
		   scalar_type & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    callCount_++;

    // compute r from system object
    system.residual(state, r_);
    // apply M
    _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), r_, 
        Mr_, is_irwls, callCount_);

    residualNorm = this->_computeNorm();

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }

    if (recomputeSystemJacobian){
      system.jacobian(state, J_);
      _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), J_, 
          MJ_, is_irwls, callCount_);
      _computeHessian();
    }

    _computeGradient();
  }

  template<typename SystemType, typename StateType>
  mpl::enable_if_t<
    ::pressio::nonlinearsolvers::constraints::system_fused_residual_jacobian<SystemType>::value
    >
  computeOperators(const SystemType & system,
		   const StateType & state,
		   scalar_type & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    callCount_++;

    system.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), r_, 
        Mr_, is_irwls, callCount_);
    residualNorm = this->_computeNorm();

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }

    if (recomputeSystemJacobian){
      _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), J_, 
          MJ_, is_irwls, callCount_);
      _computeHessian();
    }

    _computeGradient();
  }

  template< typename SystemType, typename StateType>
  mpl::enable_if_t<
    ::pressio::nonlinearsolvers::constraints::system_residual_jacobian<SystemType>::value
    >
  residualNorm(const SystemType & system,
	       const StateType & state,
	       scalar_type & residualNorm) const
  {
    system.residual(state, r_);
    _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), r_, Mr_);
    residualNorm = this->_computeNorm();

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

  template< typename SystemType, typename StateType>
  mpl::enable_if_t<
    ::pressio::nonlinearsolvers::constraints::system_fused_residual_jacobian<SystemType>::value
    >
  residualNorm(const SystemType & system,
	       const StateType & state,
	       scalar_type & residualNorm) const
  {
    // here we query system to recompute r_ only (that is why we pass false)
    system.residualAndJacobian(state, r_, J_, false);
    _applyWeightingHelper<ResidualType,JacobianType>(functorM_.get(), r_, Mr_);
    residualNorm = this->_computeNorm();

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

private:
  scalar_type _computeNorm() const
  {
    return std::sqrt(::pressio::ops::dot(r_, Mr_));
  }

  void _computeHessian()
  {
    constexpr auto beta  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::product(pT, pnT, alpha, J_, MJ_, beta, H_);
  }


  void _computeGradient()
  {
    constexpr auto beta  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto alpha = ::pressio::utils::Constants<scalar_type>::one();
    // compute gradient (g_ = J^T M r)
    ::pressio::ops::product(pT, alpha, J_, Mr_, beta, g_);
    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::Constants<scalar_type>::negOne());
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_WEIGHTED_RJ_API_HPP_
