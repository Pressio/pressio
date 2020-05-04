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

#ifndef ROM_LSPG_STEADY_SYSTEM_HPP_
#define ROM_LSPG_STEADY_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace steady{


template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type,
  typename enable = void
  >
class System;

template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
class System<
  app_type, lspg_state_type, lspg_residual_type,
  lspg_jacobian_type, residual_policy_type,
  jacobian_policy_type
  >
{

  const app_type & app_;
  const residual_policy_type & residualEvaluator_;
  const jacobian_policy_type & jacobianEvaluator_;
  mutable lspg_residual_type R_;
  mutable lspg_jacobian_type J_;

public:
  // these need to be public because are detected by solver
  using scalar_type = typename app_type::scalar_type;
  using state_type	= lspg_state_type;
  using residual_type	= lspg_residual_type;
  using jacobian_type	= lspg_jacobian_type;

public:
  System() = delete;
  ~System() = default;

  template <typename lspg_state_t>
  System(const app_type & appIn,
	 const residual_policy_type & resPolicyObj,
	 const jacobian_policy_type & jacPolicyObj,
	 lspg_state_t	& yROM)
    : app_(appIn),
      residualEvaluator_(resPolicyObj),
      jacobianEvaluator_(jacPolicyObj),
      R_( residualEvaluator_(yROM, app_) ),
      J_( jacobianEvaluator_(yROM, app_) ){}

public:
  void residual(const lspg_state_type & romState, lspg_residual_type & R) const{
    (this->residualEvaluator_).template operator()(romState, R, app_);
  }

  void jacobian(const lspg_state_type & romState, lspg_jacobian_type & J) const{
    (this->jacobianEvaluator_).template operator()(romState, J, app_);
  }

  lspg_residual_type residual(const lspg_state_type & romState) const{
    return (this->residualEvaluator_).template operator()(romState, app_);
  }

  lspg_jacobian_type jacobian(const lspg_state_type & romState) const{
    return (this->jacobianEvaluator_).template operator()(romState, app_);
  }

  scalar_type operator()(const state_type & romState) const
  {
    this->residual(romState, R_);
    constexpr auto two = ::pressio::utils::constants::two<scalar_type>();
    const auto norm    = pressio::ops::norm2(R_);
    return norm*norm;
  }

  void gradient( const state_type & romState, state_type & g) const
  {
    this->residual(romState, R_);
    this->jacobian(romState, J_);
    constexpr auto beta  = ::pressio::utils::constants::zero<scalar_type>();
    constexpr auto alpha = ::pressio::utils::constants::two<scalar_type>();
    ::pressio::ops::product(::pressio::transpose(), alpha, J_, R_, beta, g);
  }
};//end class

}}}} // end namespace pressio::ode::lspg::steady
#endif
