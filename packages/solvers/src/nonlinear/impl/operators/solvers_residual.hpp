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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUALOPERATORS_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUALOPERATORS_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename r_t>
class ResidualOperator
{
  using sc_t = typename ::pressio::containers::details::traits<r_t>::scalar_t;

  r_t r_;

  // auxR_ is used for residualNorm method so that we don't modify the real operator r_
  // which must be the same once computeOperators is called.
  mutable r_t auxR_;

public:
  ResidualOperator() = delete;
  ResidualOperator(ResidualOperator const &) = default;
  ResidualOperator & operator=(ResidualOperator const &) = default;
  ResidualOperator(ResidualOperator && o) = default;
  ResidualOperator & operator=(ResidualOperator && o) = default;
  ~ResidualOperator() = default;

  template <typename system_t, typename state_t>
  ResidualOperator(const system_t & system,
		   const state_t & state)
    : r_( system.createResidual() ),
      auxR_( system.createResidual() )
  {
    ::pressio::ops::set_zero(r_);
  }

public:
  void resetForNewCall()	{ /* no op */ }
  r_t & residualRef()		{ return r_; }
  const r_t & residualCRef() const { return r_; }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("Residual operator does not have parameters");
  }

  sc_t getParameter(std::string key) const {
    throw std::runtime_error("Residual operator does not have parameters");
    return {};
  }

  template<typename system_t, typename state_t>
  void computeOperators(const system_t & sys,
			const state_t & state,
			sc_t & residualNorm)
  {
    sys.residual(state, r_);
    residualNorm = ::pressio::ops::norm2(r_);
  }

  template< typename system_t, typename state_t>
  void residualNorm(const system_t & system,
		    const state_t & state,
		    sc_t & residualNorm) const
  {
    system.residual(state, auxR_);
    residualNorm = ::pressio::ops::norm2(auxR_);
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_RESIDUALOPERATORS_HPP_
