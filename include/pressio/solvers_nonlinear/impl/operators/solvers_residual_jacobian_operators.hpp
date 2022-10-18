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

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <typename ResidualType, typename JacobianType>
class ResidualJacobianOperators
{
public:
  using residual_norm_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));

private:
  ResidualType r_;

  // J_ is mutable because when we have the fused_res_jac api, to compute the
  // residualNorm (which is a const method), we call
  // residualAndJacobian without recomputing J_ but J_ still needs to be passed
  mutable JacobianType J_;

  // auxR_ is used for residualNorm method so that we don't modify the real operator r_
  // which must be the same once computeOperators is called.
  mutable ResidualType auxR_;

public:
  ResidualJacobianOperators() = delete;
  ResidualJacobianOperators(ResidualJacobianOperators const &) = default;
  ResidualJacobianOperators & operator=(ResidualJacobianOperators const &) = default;
  ResidualJacobianOperators(ResidualJacobianOperators && o) = default;
  ResidualJacobianOperators & operator=(ResidualJacobianOperators && o) = default;
  ~ResidualJacobianOperators() = default;

  template <
    typename SystemType/*,
    mpl::enable_if_t<
      ::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value or
      ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value,
      int
      > = 0*/
    >
  ResidualJacobianOperators(const SystemType & system)
    : r_( system.createResidual() ),
      J_( system.createJacobian() ),
      auxR_(::pressio::ops::clone(r_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
  }

public:
  void resetForNewCall()	{ /* no op */ }
  ResidualType & residualRef()		{ return r_; }
  JacobianType & jacobianRef()		{ return J_; }
  const ResidualType & residualCRef() const { return r_; }
  const JacobianType & jacobianCRef() const { return J_; }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("ResidualJacobian operators do not have parameters");
  }

  // fix the template here
  residual_norm_type getParameter(std::string key) const {
    throw std::runtime_error("ResidualJacobian operators do not have parameters");
    return {};
  }

  template<typename SystemType, typename StateType>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithResidualAndJacobian<SystemType> 
  void 
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value
  >  
#endif
  computeOperators(const SystemType & sys,
		   const StateType & state,
		   residual_norm_type & residualNorm,
		   bool recomputeSystemJacobian=true)
  {

    sys.residual(state, r_);
    residualNorm = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }

    if  (recomputeSystemJacobian){
      sys.jacobian(state, J_);
    }
  }

  template<typename SystemType, typename state_t>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithFusedResidualAndJacobian<SystemType> 
  void 
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value
  >  
#endif
  computeOperators(const SystemType & sys,
		   const state_t & state,
		   residual_norm_type & residualNorm,
		   bool recomputeSystemJacobian=true)
  {
    sys.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    residualNorm = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

  template< typename SystemType, typename state_t>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithResidualAndJacobian<SystemType> 
  void 
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithResidualAndJacobian<SystemType>::value
  >  
#endif
  residualNorm(const SystemType & system,
	       const state_t & state,
	       residual_norm_type & residualNorm) const
  {
    system.residual(state, auxR_);
    residualNorm = ::pressio::ops::norm2(auxR_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

  template< typename SystemType, typename state_t>
#ifdef PRESSIO_ENABLE_CXX20
  requires SystemWithFusedResidualAndJacobian<SystemType> 
  void 
#else
  mpl::enable_if_t<
  ::pressio::nonlinearsolvers::SystemWithFusedResidualAndJacobian<SystemType>::value
  >  
#endif
  residualNorm(const SystemType & system,
	       const state_t & state,
	       residual_norm_type & residualNorm) const
  {
    system.residualAndJacobian(state, auxR_, J_, false);
    residualNorm = ::pressio::ops::norm2(auxR_);

    if (std::isnan(residualNorm)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUAL_JACOBIAN_OPERATORS_HPP_
