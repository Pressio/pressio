/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_system.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_SYSTEM_HPP_
#define ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<
  class ScalarType,
  class LspgStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
class SteadySystem
{
  std::reference_wrapper<const ResidualPolicyType> residualPolicy_;
  std::reference_wrapper<const JacobianPolicyType> jacobianPolicy_;
  mutable LspgResidualType R_;
  mutable LspgJacobianType J_;

public:
  // these need to be public because are detected by solver
  using scalar_type	= ScalarType;
  using state_type	= LspgStateType;
  using residual_type	= LspgResidualType;
  using jacobian_type	= LspgJacobianType;

public:
  SteadySystem() = delete;
  SteadySystem(const SteadySystem &) = default;
  SteadySystem & operator=(const SteadySystem &) = delete;
  SteadySystem(SteadySystem &&) = default;
  SteadySystem & operator=(SteadySystem &&) = delete;
  ~SteadySystem() = default;

  SteadySystem(const ResidualPolicyType & resPolicyObj,
	       const JacobianPolicyType & jacPolicyObj)
    : residualPolicy_(resPolicyObj),
      jacobianPolicy_(jacPolicyObj),
      R_(residualPolicy_.get().create()),
      J_(jacobianPolicy_.get().create())
    {}

public:
  LspgResidualType createResidual() const{ return R_; }
  LspgJacobianType createJacobian() const{ return J_; }

  void residual(const LspgStateType & romState, LspgResidualType & R) const
  {
    residualPolicy_.get().template compute(romState, R);
  }

  void jacobian(const LspgStateType & romState, LspgJacobianType & J) const
  {
    jacobianPolicy_.get().template compute(romState, J);
  }

  //
  // the following are needed to interface with the optimizers
  //
  // scalar_type operator()(const state_type & romState) const
  // {
  //   // scalar_type normR = {};
  //   this->residual(romState, R_);
  //   const auto normR = ::pressio::ops::norm2(R_);
  //   return normR*normR;
  // }

  // void gradient( const state_type & romState, state_type & g) const
  // {
  //   // scalar_type normR = {};
  //   this->residual(romState, R_);//, ::pressio::Norm::L2, normR);
  //   this->jacobian(romState, J_);
  //   constexpr auto beta  = ::pressio::utils::Constants<scalar_type>::zero();
  //   constexpr auto alpha = ::pressio::utils::Constants<scalar_type>::two();
  //   ::pressio::ops::product(::pressio::transpose(), alpha, J_, R_, beta, g);
  // }
};

}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_SYSTEM_HPP_
