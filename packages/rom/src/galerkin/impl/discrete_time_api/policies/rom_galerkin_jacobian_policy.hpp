/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_jacobian_policy_discrete_time_api.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<class galerkin_jacobian_type, typename projection_policy_t>
class JacobianPolicyDiscreteTimeApi : private projection_policy_t
{

public:
  JacobianPolicyDiscreteTimeApi() = delete;
  JacobianPolicyDiscreteTimeApi(const JacobianPolicyDiscreteTimeApi &) = default;
  JacobianPolicyDiscreteTimeApi & operator=(const JacobianPolicyDiscreteTimeApi &) = delete;
  JacobianPolicyDiscreteTimeApi(JacobianPolicyDiscreteTimeApi &&) = default;
  JacobianPolicyDiscreteTimeApi & operator=(JacobianPolicyDiscreteTimeApi &&) = delete;
  ~JacobianPolicyDiscreteTimeApi() = default;

  template<typename ...Args>
  JacobianPolicyDiscreteTimeApi(std::size_t romSize,
				Args && ...args)
    : projection_policy_t(std::forward<Args>(args)...),
      romSize_(romSize){}

public:
  template <typename fom_system_t>
  galerkin_jacobian_type create(const fom_system_t & fomSystemObj) const
  {
    return galerkin_jacobian_type(romSize_, romSize_);
  }

  template <
    typename stepper_tag,
    typename galerkin_state_t,
    typename prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const galerkin_state_t & galerkinState,
	       const prev_states_t & galerkinPrevStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & timeStep,
	       galerkin_jacobian_type & galerkinJacobian) const
  {
    projection_policy_t::compute(galerkinJacobian, galerkinState,
				 fomSystemObj, galerkinPrevStates,
				 time, dt, timeStep);
  }

protected:
  const std::size_t romSize_ = {};
};

}}}}
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_JACOBIAN_POLICY_DISCRETE_TIME_API_HPP_
