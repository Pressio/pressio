/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_system.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

namespace pressio{ namespace rom{

template<
  typename app_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
class LSPGSteadySystem<
  app_type, lspg_state_type, lspg_residual_type,
  lspg_jacobian_type, residual_policy_type,
  jacobian_policy_type
  >{

  const app_type & app_;
  const residual_policy_type & residualEvaluator_;
  const jacobian_policy_type & jacobianEvaluator_;

public:
  // these need to be public because are detected by solver
  using scalar_type = typename app_type::scalar_type;
  using state_type	= lspg_state_type;
  using residual_type	= lspg_residual_type;
  using jacobian_type	= lspg_jacobian_type;

public:
  LSPGSteadySystem() = delete;
  ~LSPGSteadySystem() = default;

  LSPGSteadySystem(const app_type & appIn,
		   const residual_policy_type & resPolicyObj,
		   const jacobian_policy_type & jacPolicyObj)
    : app_(appIn),
      residualEvaluator_(resPolicyObj),
      jacobianEvaluator_(jacPolicyObj){}

public:
  void residual(const lspg_state_type & y,
		lspg_residual_type & R) const{
    (this->residualEvaluator_).template operator()(y, R, app_);
  }

  void jacobian(const lspg_state_type & y,
		lspg_jacobian_type & J) const{
    (this->jacobianEvaluator_).template operator()(y, J, app_);
  }

  lspg_residual_type residual(const lspg_state_type & y) const{
    return (this->residualEvaluator_).template operator()(y, app_);
  }

  lspg_jacobian_type jacobian(const lspg_state_type & y) const{
    return (this->jacobianEvaluator_).template operator()(y, app_);
  }
};//end class

}} // end namespace pressio::ode
#endif
