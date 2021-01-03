/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_bdf2_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_BDF2_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_BDF2_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{

template<
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename system_type,
  typename aux_stepper_t,
  typename residual_policy_t,
  typename jacobian_policy_t,
  bool policies_are_standard
  >
class StepperBDF2
{

public:
  // these need to be here because are detected by solver
  using scalar_type = scalar_t;
  using state_type  = state_t;
  using residual_type = residual_t;
  using jacobian_type = jacobian_t;

  using tag_name = ::pressio::ode::implicitmethods::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr types::stepper_order_t order_value = 2;
  static constexpr std::size_t numAuxStates = 2;
  using stencil_states_t = StencilStatesManager<state_t, numAuxStates>;

private:
  // auxiliary stepper
  aux_stepper_t & auxStepper_;

  scalar_t rhsEvaluationTime_  = {};
  scalar_t dt_ = {};
  types::step_t stepNumber_  = {};
  std::reference_wrapper<const system_type> systemObj_;
  stencil_states_t stencilStates_;

  // state object to ensure the strong guarantee for handling excpetions
  state_type recoveryState_;

  // policies
  ::pressio::utils::instance_or_reference_wrapper<residual_policy_t> resPolicy_;
  ::pressio::utils::instance_or_reference_wrapper<jacobian_policy_t> jacPolicy_;

public:
  StepperBDF2() = delete;
  StepperBDF2(const StepperBDF2 & other)  = default;
  StepperBDF2 & operator=(const StepperBDF2 & other) = delete;
  StepperBDF2(StepperBDF2 && other)  = default;
  StepperBDF2 & operator=(StepperBDF2 && other) = delete;
  ~StepperBDF2() = default;

  StepperBDF2(const state_type & state,
	      const system_type & systemObj,
	      const mpl::remove_cvref_t<residual_policy_t> & resPolicyObj,
	      const mpl::remove_cvref_t<jacobian_policy_t> & jacPolicyObj,
	      aux_stepper_t & auxStepper)
    : auxStepper_{auxStepper},
      systemObj_{systemObj},
      stencilStates_{state},
      recoveryState_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
  {}

  template <
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int> = 0
    >
  StepperBDF2(const state_type & state,
	      const system_type & systemObj,
	      aux_stepper_t & auxStepper)
    : auxStepper_{auxStepper},
      systemObj_{systemObj},
      stencilStates_{state},
      recoveryState_{state},
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
    PRESSIOLOG_DEBUG("bdf2 stepper: do step");

    auto dummyGuesser =
      [](const types::step_t &, const scalar_t &, state_type &)
      { /*no op*/ };

    doStepImpl(odeState, currentTime, dt, stepNumber, solver, dummyGuesser);
  }

  template<typename solver_type, typename guess_callback_t>
  void doStep(state_type & odeState,
	      const scalar_t &  currentTime,
	      const scalar_t &  dt,
	      const types::step_t & stepNumber,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    PRESSIOLOG_DEBUG("bdf2 stepper: do step with callback to state guesser");
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
      <tag_name>(odeState,  stencilStates_,
		 systemObj_.get(),
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
  void doStepImpl(state_type & odeState,
		  const scalar_t & currentTime,
		  const scalar_t & dt,
		  const types::step_t & stepNumber,
		  solver_type & solver,
		  guess_callback_t && guesserCb)
  {
    static_assert(::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF2 stepper");

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeSolution.

      for bdf1, this is done by solving:  y_n+1 = y_n + dt*f(t_n+1, y_n+1)
     */

    dt_ = dt;
    rhsEvaluationTime_ = currentTime + dt;
    stepNumber_ = stepNumber;

    // first step, use auxiliary stepper
    if (stepNumber == 1){
      // step ==1 means that we are going from y_0 to y_1
      // stencilStates_(0) now holds y_0
      ::pressio::ops::deep_copy(stencilStates_.stateAt(ode::n()), odeState);
      ::pressio::ops::deep_copy(stencilStates_.stateAt(ode::nMinusOne()), odeState);
      auxStepper_.doStep(odeState, currentTime, dt, stepNumber, solver);
    }
    if (stepNumber >= 2)
    {
      /*
	at step == 2 we are going from t_1 to t_2, so
	 we want to compute:
	 y_2 = (4/3)*y_1 - (1/3)*y_0 + (2/3)*dt*f(y_2, t_2)

	 upon entering step=2, we have:
	 odeState	contains y_1 computed from auxiliary stepper
	 auxStates(n()) contains y_0

	 so we need to update data such that
	 n becomes nm1 and odeState becomes n

	step == 3 means that we are going from t_2 to t_3, so
	we want to compute:
	y_3 = (4/3)*y_2 - (1/3)*y_1 + (2/3)*dt*f(y_3, t_3)
	so this means that: y_n = y_2 and y_n-1 = y_1

	and so on...
      */

      auto & odeState_n   = stencilStates_.stateAt(ode::n());
      auto & odeState_nm1 = stencilStates_.stateAt(ode::nMinusOne());
      ::pressio::ops::deep_copy(recoveryState_, odeState_nm1);
      ::pressio::ops::deep_copy(odeState_nm1, odeState_n);
      ::pressio::ops::deep_copy(odeState_n, odeState);

      try{
	guesserCb(stepNumber, rhsEvaluationTime_, odeState);
	solver.solve(*this, odeState);
      }
      catch (::pressio::eh::nonlinear_solve_failure const & e)
	{
	  ::pressio::ops::deep_copy(odeState, odeState_n);
	  ::pressio::ops::deep_copy(odeState_n, odeState_nm1);
	  ::pressio::ops::deep_copy(odeState_nm1, recoveryState_);

	  // now throw
	  throw ::pressio::eh::time_step_failure();
	}
    }
  }
};

}}} // end namespace pressio::ode::implicitmethods
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_BDF2_IMPL_HPP_
