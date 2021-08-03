/*
//@HEADER
// ************************************************************************
//
// solvers_weighting_irwls.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_WEIGHTING_IRWLS_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_WEIGHTING_IRWLS_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template <class ResidualType, class JacobType, class ScalarType>
class IrwWeightingOperator
{
public:
  using scalar_type = ScalarType;

private:
  mutable ResidualType w_;
  scalar_type p_ = ::pressio::utils::Constants<scalar_type>::one();
  scalar_type exponent_ = {};

public:
  IrwWeightingOperator() = delete;
  IrwWeightingOperator(IrwWeightingOperator const &) = default;
  IrwWeightingOperator & operator=(IrwWeightingOperator const &) = default;
  IrwWeightingOperator(IrwWeightingOperator && o) = default;
  IrwWeightingOperator & operator=(IrwWeightingOperator && o) = default;
  ~IrwWeightingOperator() = default;

  template <
    typename SystemType,
    mpl::enable_if_t<
    ::pressio::nonlinearsolvers::constraints::system_residual_jacobian<SystemType>::value or
    ::pressio::nonlinearsolvers::constraints::system_fused_residual_jacobian<SystemType>::value,
    int > = 0
    >
  IrwWeightingOperator(const SystemType & system)
    : w_(system.createResidual())
  {
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::fill(w_, one);
    this->computeExponent();
  }

public:
  void set_p(scalar_type pIn)
  {
    p_ = pIn;
    this->computeExponent();
  }

  void operator()(const ResidualType & rIn, ResidualType & Wr) const
  {
    this->computeWeights(rIn);
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::elementwise_multiply(one, w_, rIn, zero, Wr);
  }

  void operator()(const JacobType & Jin, JacobType & WJ) const
  {
    // don't compute weights here since they have been computed above

    // view the weights as a diagonal matrix
    const auto wMat = ::pressio::as_diagonal_matrix(w_);

    // WJ = W * Jin
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::product(::pressio::nontranspose(),
			    ::pressio::nontranspose(),
			    one, wMat, Jin,
			    zero, WJ);
  }

private:
  void computeExponent()
  {
    constexpr auto two = ::pressio::utils::Constants<scalar_type>::two();
    if (p_==two) throw std::runtime_error("irwls does not support using p=2!");
    exponent_ = (p_ - two);
  }

  void computeWeights(const ResidualType & err) const
  {
    // compute w = |e|^(p-2)
    // use a small epsilong to guard against diving by zero
    // when exponent < 0

    constexpr auto two = ::pressio::utils::Constants<scalar_type>::two();
    if (p_ > two){
      ::pressio::ops::abs_pow(w_, err, exponent_);
    }
    else{
      ::pressio::ops::abs_pow(w_, err, exponent_, 0.00001);
    }
  }
};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_WEIGHTING_IRWLS_HPP_
