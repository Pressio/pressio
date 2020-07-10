/*
//@HEADER
// ************************************************************************
//
// rom_preconditioned_decorator_jacobian.hpp
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

#ifndef ROM_PRECONDITIONED_DECORATOR_JACOBIAN_HPP_
#define ROM_PRECONDITIONED_DECORATOR_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace decorator{

template <typename preconditionable_policy>
class PreconditionedJacobianPolicy : public preconditionable_policy
{

public:
  using typename preconditionable_policy::apply_jac_return_t;
  using preconditionable_policy::fomStatesMngr_;

public:
  PreconditionedJacobianPolicy() = delete;
  PreconditionedJacobianPolicy(const preconditionable_policy & obj) 
    : preconditionable_policy(obj)
  {}

  template <typename ... Args>
  PreconditionedJacobianPolicy(Args && ... args)
    : preconditionable_policy(std::forward<Args>(args)...)
  {}

  ~PreconditionedJacobianPolicy() = default;

public:
  template <typename app_t>
  apply_jac_return_t create(const app_t & app) const
  {
    return preconditionable_policy::create(app);
  }

  //-------------------------------
  // unsteady
  //-------------------------------
  template <
    typename stepper_tag,
    typename prev_states_mgr,
    typename ode_state_t,
    typename app_t,
    typename scalar_t,
    typename ode_jac_t
  >
  void compute(const ode_state_t & odeState,
      const prev_states_mgr & prevStatesMgr,
  		const app_t & app,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const ::pressio::ode::types::step_t & step,
		  ode_jac_t & odeJacobian) const
  {
    preconditionable_policy::template compute<stepper_tag>(odeState, prevStatesMgr, app, t, dt, step, odeJacobian);
    const auto & yFom = fomStatesMngr_.getCRefToCurrentFomState();
    app.applyPreconditioner(*yFom.data(), t, *odeJacobian.data());
  }

  //-------------------------------
  // steady 
  //-------------------------------
  template <typename state_t, typename jac_t, typename app_t>
  void compute(const state_t & stateObj,
               jac_t & jacob,
               const app_t & app) const
  {
    preconditionable_policy::compute(stateObj, jacob, app);
    const auto & yFom = fomStatesMngr_.getCRefToCurrentFomState();
    app.applyPreconditioner(*yFom.data(), *jacob.data());
  }

};//end class

}}} //end namespace pressio::rom::decorator

#endif
