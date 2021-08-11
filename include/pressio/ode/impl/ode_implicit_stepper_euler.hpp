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

namespace pressio{ namespace ode{ namespace impl{

template<
  typename ScalarType,
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename ResidualPolicyType,
  typename JacobianPolicyType,
  bool policies_are_standard
  >
class StepperBDF1
{

public:
  // these aliases are detected by solver
  using scalar_type	= ScalarType;
  using state_type	= StateType;
  using residual_type	= ResidualType;
  using jacobian_type	= JacobianType;

  using tag_name = ::pressio::ode::implicitmethods::BDF1;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr stepper_order_type order_value = 1;

private:
  ScalarType rhsEvaluationTime_  = {};
  ScalarType dt_ = {};
  ::pressio::ode::step_count_type stepNumber_  = {};
  std::reference_wrapper<const SystemType> systemObj_;
  ImplicitStencilStatesContainer<StateType, 1> stencilStates_;
  ::pressio::utils::InstanceOrReferenceWrapper<ResidualPolicyType> resPolicy_;
  ::pressio::utils::InstanceOrReferenceWrapper<JacobianPolicyType> jacPolicy_;

public:
  StepperBDF1() = delete;
  StepperBDF1(const StepperBDF1 & other)  = default;
  StepperBDF1 & operator=(const StepperBDF1 & other) = delete;
  StepperBDF1(StepperBDF1 && other)  = default;
  StepperBDF1 & operator=(StepperBDF1 && other) = delete;
  ~StepperBDF1() = default;

  template<
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int > = 0
    >
  StepperBDF1(const state_type & state,
	      const SystemType & systemObj)
    : systemObj_{systemObj},
      stencilStates_(state), //stencilstates handles right semantics
      resPolicy_{},
      jacPolicy_{}
  {}

  StepperBDF1(const state_type & state,
	      const SystemType & systemObj,
	      ResidualPolicyType && resPolicyObj,
	      JacobianPolicyType && jacPolicyObj)
    : systemObj_{systemObj},
      stencilStates_(state), //stencilstates handles right semantics
      resPolicy_{std::forward<ResidualPolicyType>(resPolicyObj)},
      jacPolicy_{std::forward<JacobianPolicyType>(jacPolicyObj)}
  {}

  template<
    class _ResidualPolicyType = ResidualPolicyType,
    class _JacobianPolicyType = JacobianPolicyType,
    mpl::enable_if_t<!std::is_reference<_ResidualPolicyType>::value and
		     !std::is_reference<_JacobianPolicyType>::value, int> = 0
    >
  StepperBDF1(const state_type & state,
	      const SystemType & systemObj,
	      const _ResidualPolicyType & resPolicyObj,
	      const _JacobianPolicyType & jacPolicyObj)
    : systemObj_{systemObj},
      stencilStates_(state), //stencilstates handles right semantics
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
  {}

public:
  ::pressio::ode::stepper_order_type order() const
  {
    return order_value;
  }

  template<typename solver_type, typename ...Args>
  void doStep(state_type & odeState,
	      const ScalarType & currentTime,
	      const ScalarType & dt,
	      const ::pressio::ode::step_count_type & stepNumber,
	      solver_type & solver,
	      Args&& ...args)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step");

    auto dummyGuesser =
      [](const ::pressio::ode::step_count_type &,
         const ScalarType &,
         state_type &)
      { /*no op*/ };

    doStepImpl(odeState, currentTime, dt, stepNumber,
	       dummyGuesser, solver,
	       std::forward<Args>(args)...);
  }

  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void doStep(state_type & odeState,
	      const ScalarType & currentTime,
	      const ScalarType & dt,
	      const ::pressio::ode::step_count_type & stepNumber,
	      guess_callback_t && guesserCb,
	      solver_type & solver,
	      Args&& ...args)
  {
    PRESSIOLOG_DEBUG("bdf1 stepper: do step with callback to state guesser");
    doStepImpl(odeState, currentTime, dt, stepNumber,
	       guesserCb, solver, std::forward<Args>(args)...);
  }

  ResidualType createResidual() const
  {
    return resPolicy_.get().create(systemObj_.get());
  }

  JacobianType createJacobian() const
  {
    return jacPolicy_.get().create(systemObj_.get());
  }

  void residual(const StateType & odeState, ResidualType & R) const
  {
    resPolicy_.get().template compute
      <tag_name>(odeState, stencilStates_, systemObj_.get(),
		 rhsEvaluationTime_, dt_, stepNumber_, R);
  }

  void jacobian(const StateType & odeState, JacobianType & J) const
  {
    jacPolicy_.get().template compute<
      tag_name>(odeState, stencilStates_, systemObj_.get(),
                rhsEvaluationTime_, dt_, stepNumber_, J);
  }

private:
  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void doStepImpl(state_type & odeSolution,
		  const ScalarType & currentTime,
		  const ScalarType & dt,
		  const ::pressio::ode::step_count_type & stepNumber,
		  guess_callback_t && guesserCb,
		  solver_type & solver,
		  Args&& ...args)
  {
    static_assert(::pressio::ode::legitimate_solver_for_implicit_stepper<
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
    auto & odeState_n = stencilStates_(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeSolution);

    // if provided, callback to provide a guess for the odeSolution
    guesserCb(stepNumber, rhsEvaluationTime_, odeSolution);

    try{
      solver.solve(*this, odeSolution, std::forward<Args>(args)...);
    }
    catch (::pressio::eh::NonlinearSolveFailure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeSolution to that
      auto & rollBackState = stencilStates_(ode::n());
      ::pressio::ops::deep_copy(odeSolution, rollBackState);

      // now throw
      throw ::pressio::eh::TimeStepFailure();
    }
  }
};

}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_EULER_IMPL_HPP_
