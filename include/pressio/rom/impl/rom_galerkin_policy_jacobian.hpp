/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_policy_jacobian.hpp
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

#ifndef ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<class GalerkinJacobianType, typename ProjectionPolicyType>
class JacobianPolicy : private ProjectionPolicyType
{

public:
  // required
  using jacobian_type = GalerkinJacobianType;

public:
  JacobianPolicy() = delete;
  JacobianPolicy(const JacobianPolicy &) = default;
  JacobianPolicy & operator=(const JacobianPolicy &) = delete;
  JacobianPolicy(JacobianPolicy &&) = default;
  JacobianPolicy & operator=(JacobianPolicy &&) = delete;
  ~JacobianPolicy() = default;

  template<typename ...Args>
  JacobianPolicy(Args && ...args)
    : ProjectionPolicyType(std::forward<Args>(args)...){}

public:
  jacobian_type create() const
  {
    return ProjectionPolicyType::template create<jacobian_type>();
  }

  template <
    class GalerkinStateType,
    class StencilStatesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::StepScheme name,
		  const GalerkinStateType & galerkinState,
		  const StencilStatesContainerType & stencilStates,
		  const ScalarType & time_np1,
		  const ScalarType & dt,
		  const StepType & currentStepNumber,
		  jacobian_type & galerkinJacobian) const
  {
    ProjectionPolicyType::compute(galerkinJacobian, galerkinState, time_np1);

    if (name == ::pressio::ode::StepScheme::BDF1){
      ::pressio::ode::impl::discrete_time_jacobian(galerkinJacobian, dt,
						   ode::BDF1());
    }
    else if (name == ::pressio::ode::StepScheme::BDF2){
      ::pressio::ode::impl::discrete_time_jacobian(galerkinJacobian, dt,
						   ode::BDF2());
    }
    else if (name == ::pressio::ode::StepScheme::CrankNicolson){
      ::pressio::ode::impl::discrete_time_jacobian(galerkinJacobian, dt,
						   ode::CrankNicolson());
    }
  }

};

}}}}
#endif  // ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_HPP_
