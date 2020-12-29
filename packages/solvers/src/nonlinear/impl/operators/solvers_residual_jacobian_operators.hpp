/*
//@HEADER
// ************************************************************************
//
// solvers_residual_jacobian_operators.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename r_t, typename j_t>
class ResidualJacobianOperators
{
  using sc_t = typename r_t::traits::scalar_t;

  r_t r_;

  // J_ is mutable because when we have the fused_res_jac api, to compute the
  // residualNorm (which is a const method), we call
  // residualAndJacobian without recomputing J_ but J_ still needs to be passed
  mutable j_t J_;

  // auxR_ is used for residualNorm method so that we don't modify the real operator r_
  // which must be the same once computeOperators is called.
  mutable r_t auxR_;

public:
  ResidualJacobianOperators() = delete;
  ResidualJacobianOperators(ResidualJacobianOperators const &) = default;
  ResidualJacobianOperators & operator=(ResidualJacobianOperators const &) = default;
  ResidualJacobianOperators(ResidualJacobianOperators && o) = default;
  ResidualJacobianOperators & operator=(ResidualJacobianOperators && o) = default;
  ~ResidualJacobianOperators() = default;

  template <
    typename system_t, typename state_t,
    mpl::enable_if_t<
      pressio::solvers::constraints::system_residual_jacobian<system_t>::value or
      pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value,
      int
     > = 0
  >
  ResidualJacobianOperators(const system_t & system,
			    const state_t & state)
    : r_( system.createResidual() ),
      J_( system.createJacobian() ),
      auxR_( system.createResidual() )
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
  }

public:
  void resetForNewCall()	{ /* no op */ }
  r_t & residualRef()		{ return r_; }
  j_t & jacobianRef()		{ return J_; }
  const r_t & residualCRef() const { return r_; }
  const j_t & jacobianCRef() const { return J_; }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("ResidualJacobian operators do not have parameters");
  }

  sc_t getParameter(std::string key) const {
    throw std::runtime_error("ResidualJacobian operators do not have parameters");
    return {};
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_residual_jacobian<system_t>::value
  >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian=true)
  {

    sys.residual(state, r_);
    residualNorm = ::pressio::ops::norm2(r_);

    if  (recomputeSystemJacobian){
      sys.jacobian(state, J_);
    }
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian=true)
  {
    sys.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    residualNorm = ::pressio::ops::norm2(r_);
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system,
         const state_t & state,
         sc_t & residualNorm) const
  {
    system.residual(state, auxR_);
    residualNorm = ::pressio::ops::norm2(auxR_);
  }

  template< typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_residual_jacobian<system_t>::value
    >
  residualNorm(const system_t & system,
         const state_t & state,
         sc_t & residualNorm) const
  {
    system.residualAndJacobian(state, auxR_, J_, false);
    residualNorm = ::pressio::ops::norm2(auxR_);
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
