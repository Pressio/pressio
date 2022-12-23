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

template <class ResidualType, class JacobianType>
class FusedResidualJacobianOperators
{
public:
  using residual_norm_value_type = decltype(::pressio::ops::norm2(std::declval<ResidualType>()));

private:
  ResidualType r_;
  mutable JacobianType J_;
  mutable ResidualType auxR_;
  residual_norm_value_type residualNorm_;

public:
  template <class SystemType>
  FusedResidualJacobianOperators(const SystemType & system)
    : r_( system.createResidual() ),
      J_( system.createJacobian() ),
      auxR_(::pressio::ops::clone(r_))
  {
    ::pressio::ops::set_zero(r_);
    ::pressio::ops::set_zero(J_);
  }

public:
  const ResidualType & residualCRef() const { return r_; }
  const JacobianType & jacobianCRef() const { return J_; }

  template<class SystemType, class state_t>
  void compute(const SystemType & sys,
	       const state_t & state,
	       bool recomputeSystemJacobian=true)
  {
    sys.residualAndJacobian(state, r_, J_, recomputeSystemJacobian);
    residualNorm_ = ::pressio::ops::norm2(r_);

    if (std::isnan(residualNorm_)){
      throw ::pressio::eh::ResidualHasNans();
    }
  }

  residual_norm_value_type residualNorm() const{
    return residualNorm_:
  }

  template< class SystemType, class state_t>
  void residualNorm(const SystemType & system,
		    const state_t & state,
		    residual_norm_value_type & residualNorm) const
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
