/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_euler_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename system_type,
  typename residual_policy_t,
  typename jacobian_policy_t,
  bool policies_are_standard
  >
class StepperBDF1
{

public:
  // these aliases are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= state_t;
  using residual_type	= residual_t;
  using jacobian_type	= jacobian_t;

  using tag_name = ::pressio::ode::implicitmethods::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr types::stepper_order_t order_value = 1;
  static constexpr std::size_t numAuxStates = 1;
  using stencil_states_t = StencilStatesManager<state_t, numAuxStates>;

private:
  scalar_t rhsEvaluationTime_  = {};
  scalar_t dt_ = {};
  types::step_t stepNumber_  = {};
  std::reference_wrapper<const system_type> systemObj_;
  stencil_states_t stencilStates_;
  ::pressio::utils::instance_or_reference_wrapper<residual_policy_t> resPolicy_;
  ::pressio::utils::instance_or_reference_wrapper<jacobian_policy_t> jacPolicy_;

public:
  StepperBDF1() = delete;
  StepperBDF1(const StepperBDF1 & other)  = default;
  StepperBDF1 & operator=(const StepperBDF1 & other) = delete;
  StepperBDF1(StepperBDF1 && other)  = default;
  StepperBDF1 & operator=(StepperBDF1 && other) = delete;
  ~StepperBDF1() = default;

  StepperBDF1(const state_type & state,
	      const system_type & systemObj,
	      const mpl::remove_cvref_t<residual_policy_t> & resPolicyObj,
	      const mpl::remove_cvref_t<jacobian_policy_t> & jacPolicyObj)
    : systemObj_{systemObj},
      stencilStates_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
  {}

  template<
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int > = 0
    >
  StepperBDF1(const state_type & state,
	      const system_type & systemObj)
    : systemObj_{systemObj},
      stencilStates_{state},
      resPolicy_{},
      jacPolicy_{}
  {}

public:
  ::pressio::ode::types::stepper_order_t order() const
  {
    return order_value;
  }

  template<typename solver_type>
  void doStep(state_type & odeState,
	      const scalar_t & currentTime,
	      const scalar_t & dt,
	      const types::step_t & stepNumber,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step");

    auto dummyGuesser =
      [](const types::step_t &, const scalar_t &, state_type &)
      { /*no op*/ };

    doStepImpl(odeState, currentTime, dt, stepNumber, solver, dummyGuesser);
  }

  template<typename solver_type, typename guess_callback_t>
  void doStep(state_type & odeState,
	      const scalar_t & currentTime,
	      const scalar_t & dt,
	      const types::step_t & stepNumber,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step with callback to state guesser");
    doStepImpl(odeState, currentTime, dt, stepNumber, solver, guesserCb);
  }

  residual_t createResidual() const
  {
    return resPolicy_.get().create(systemObj_.get());
  }

  jacobian_t createJacobian() const
  {
    return jacPolicy_.get().create(systemObj_.get());
  }

  void residual(const state_t & odeState, residual_t & R) const
  {
    resPolicy_.get().template compute
      <tag_name>(odeState, stencilStates_, systemObj_.get(),
		 rhsEvaluationTime_, dt_, stepNumber_, R);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const
  {
    jacPolicy_.get().template compute<
      tag_name>(odeState, stencilStates_, systemObj_.get(),
                rhsEvaluationTime_, dt_, stepNumber_, J);
  }

private:
  template<typename solver_type, typename guess_callback_t>
  void doStepImpl(state_type & odeSolution,
		  const scalar_t & currentTime,
		  const scalar_t & dt,
		  const types::step_t & stepNumber,
		  solver_type & solver,
		  guess_callback_t && guesserCb)
  {
    static_assert(::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF1 stepper");

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeSolution.

      for bdf1, this is done by solving:  y_n+1 = y_n + dt*f(t_n+1, y_n+1)
     */

    dt_ = dt;
    rhsEvaluationTime_ = currentTime + dt_;
    stepNumber_ = stepNumber;

    // copy current solution into y_n
    auto & odeState_n = stencilStates_.stateAt(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeSolution);

    // if provided, callback to provide a guess for the odeSolution
    guesserCb(stepNumber, rhsEvaluationTime_, odeSolution);

    try{
      solver.solve(*this, odeSolution);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeSolution to that
      auto & rollBackState = stencilStates_.stateAt(ode::n());
      ::pressio::ops::deep_copy(odeSolution, rollBackState);

      // now throw
      throw ::pressio::eh::time_step_failure();
    }
  }
};

}}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_
