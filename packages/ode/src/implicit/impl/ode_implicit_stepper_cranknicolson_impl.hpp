/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_cranknicolson_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_CRANKNICOLSON_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_CRANKNICOLSON_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{

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
class StepperCrankNicolson
{

public:
  // these need to be here because are detected by solver
  using scalar_type = scalar_t;
  using state_type  = state_t;
  using residual_type = residual_t;
  using jacobian_type = jacobian_t;

  using tag_name = ::pressio::ode::implicitmethods::CrankNicolson;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr types::stepper_order_t order_value = 2;
  using stencil_states_t = StencilStatesManager<state_t, 1>;
  using stencil_velocities_t = StencilVelocitiesManager<residual_t, 2>;

private:
  scalar_t t_np1_  = {};
  scalar_t dt_ = {};
  types::step_t stepNumber_  = {};
  std::reference_wrapper<const system_type> systemObj_;

  // state object to ensure the strong guarantee for handling excpetions
  state_type recoveryState_;

  // policies
  ::pressio::utils::instance_or_reference_wrapper<residual_policy_t> resPolicy_;
  ::pressio::utils::instance_or_reference_wrapper<jacobian_policy_t> jacPolicy_;

  // stencilStates contains: y_n
  stencil_states_t stencilStates_;
  // stencilVelocities contains f(y_n,t_n) and f(y_np1, t_np1)
  mutable stencil_velocities_t stencilVelocities_;

public:
  StepperCrankNicolson() = delete;
  StepperCrankNicolson(const StepperCrankNicolson & other)  = default;
  StepperCrankNicolson & operator=(const StepperCrankNicolson & other) = delete;
  StepperCrankNicolson(StepperCrankNicolson && other)  = default;
  StepperCrankNicolson & operator=(StepperCrankNicolson && other) = delete;
  ~StepperCrankNicolson() = default;

  StepperCrankNicolson(const state_type & state,
		       const system_type & systemObj,
		       const mpl::remove_cvref_t<residual_policy_t> & resPolicyObj,
		       const mpl::remove_cvref_t<jacobian_policy_t> & jacPolicyObj)
    : systemObj_{systemObj},
      recoveryState_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj},
      stencilStates_(state),
      stencilVelocities_(resPolicy_.get().create(systemObj))
  {}

  template <
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int > = 0
    >
  StepperCrankNicolson(const state_type & state,
		       const system_type & systemObj)
    : systemObj_{systemObj},
      recoveryState_{state},
      resPolicy_{},
      jacPolicy_{},
      stencilStates_(state),
      stencilVelocities_(resPolicy_.get().create(systemObj))
  {}

public:
  ::pressio::ode::types::stepper_order_t order() const{
    return order_value;
  }

  template<typename solver_type>
  void doStep(state_type & odeState,
	      const scalar_t &  currentTime,
	      const scalar_t &  dt,
	      const types::step_t & stepNumber,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("crankNicolson stepper: do step");

    auto dummyGuesser =
      [](const types::step_t &, const scalar_t &, state_type &)
      { /*no op*/ };

    doStepImpl(odeState, currentTime, dt, stepNumber, solver, dummyGuesser);
  }

  // overload for when we have a guesser callback
  template<typename solver_type, typename guess_callback_t>
  void doStep(state_type & odeState,
	      const scalar_t &  currentTime,
	      const scalar_t &  dt,
	      const types::step_t & stepNumber,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    PRESSIOLOG_DEBUG("crankNicolson stepper: do step with callback to state guesser");
    doStepImpl(odeState, currentTime, dt, stepNumber, solver, guesserCb);
  }

  residual_t createResidual() const{
    return resPolicy_.get().create(systemObj_.get());
  }

  jacobian_t createJacobian() const{
    return jacPolicy_.get().create(systemObj_.get());
  }

  void residual(const state_t & odeState, residual_t & R) const{
    resPolicy_.get().template compute
      <tag_name>(odeState, stencilStates_,
		 systemObj_.get(), t_np1_,
		 dt_, stepNumber_,
		 stencilVelocities_, R);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const{
    jacPolicy_.get().template compute<
      tag_name>(odeState, stencilStates_,
		systemObj_.get(), t_np1_,
		dt_, stepNumber_, J);
  }

private:
  template<typename solver_type, typename guess_callback_t>
  void doStepImpl(state_type & odeSolution,
		  const scalar_t &  currentTime,
		  const scalar_t &  dt,
		  const types::step_t & stepNumber,
		  solver_type & solver,
		  guess_callback_t && guesserCb)
  {
    static_assert(::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for CrankNicolson stepper");

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeSolution.

      for CN, this is done by solving:
      y_n+1 = y_n + 0.5*dt*[ f(t_n+1, y_n+1) + f(t_n, y_n) ]
    */

    dt_ = dt;
    t_np1_ = currentTime + dt_;
    stepNumber_ = stepNumber;

    // current solution becomes y_n
    auto & odeState_n = stencilStates_.stateAt(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeSolution);

    // if provided, callback to provide a guess for the odeSolution
    guesserCb(stepNumber, t_np1_, odeSolution);

    try{
      solver.solve(*this, odeSolution);

      // auto & fn   = stencilStates_.rhsAt(ode::n());
      // auto & fnp1 = stencilStates_.rhsAt(ode::nPlusOne());
      // //::pressio::ops::deep_copy(fn, fnp1);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      auto & rollBackState = stencilStates_.stateAt(ode::n());
      ::pressio::ops::deep_copy(odeSolution, rollBackState);
      throw ::pressio::eh::time_step_failure();
    }
  }
};

}}} // end namespace pressio::ode::implicitmethods
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_CRANKNICOLSON_IMPL_HPP_
