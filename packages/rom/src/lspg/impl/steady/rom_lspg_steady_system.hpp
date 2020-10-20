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

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template<
  typename scalar_t,
  typename fom_system_type,
  typename lspg_state_type,
  typename lspg_residual_type,
  typename lspg_jacobian_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
class System
{

  std::reference_wrapper<const fom_system_type>	     fomSystemObj_;
  std::reference_wrapper<const residual_policy_type> residualEvaluator_;
  std::reference_wrapper<const jacobian_policy_type> jacobianEvaluator_;
  mutable lspg_residual_type R_;
  mutable lspg_jacobian_type J_;

public:
  // these need to be public because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= lspg_state_type;
  using residual_type	= lspg_residual_type;
  using jacobian_type	= lspg_jacobian_type;

public:
  System() = delete;
  System(const System &) = default;
  System & operator=(const System &) = delete;
  System(System &&) = default;
  System & operator=(System &&) = delete;
  ~System() = default;

  System(const fom_system_type & fomSystemObj,
	 const residual_policy_type & resPolicyObj,
	 const jacobian_policy_type & jacPolicyObj)
    : fomSystemObj_(fomSystemObj),
      residualEvaluator_(resPolicyObj),
      jacobianEvaluator_(jacPolicyObj),
      R_(residualEvaluator_.get().create(fomSystemObj_.get())),
      J_(jacobianEvaluator_.get().create(fomSystemObj_.get()))
    {}

public:
  lspg_residual_type createResidual() const
  {
    return R_;
  }

  lspg_jacobian_type createJacobian() const
  {
    return J_;
  }

  void residual(const lspg_state_type & romState,
		lspg_residual_type & R) const
  {
    residualEvaluator_.get().template compute(romState, R, fomSystemObj_.get());
  }

  void jacobian(const lspg_state_type & romState, lspg_jacobian_type & J) const
  {
    jacobianEvaluator_.get().template compute(romState, J, fomSystemObj_.get());
  }

  // the following are needed to interface with the optimizers
  scalar_type operator()(const state_type & romState) const
  {
    // scalar_type normR = {};
    this->residual(romState, R_);
    const auto normR = ::pressio::ops::norm2(R_);
    return normR*normR;
  }

  void gradient( const state_type & romState, state_type & g) const
  {
    // scalar_type normR = {};
    this->residual(romState, R_);//, ::pressio::Norm::L2, normR);
    this->jacobian(romState, J_);
    constexpr auto beta  = ::pressio::utils::constants<scalar_type>::zero();
    constexpr auto alpha = ::pressio::utils::constants<scalar_type>::two();
    ::pressio::ops::product(::pressio::transpose(), alpha, J_, R_, beta, g);
  }
};//end class

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_LSPG_STEADY_SYSTEM_HPP_
