/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_standard.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_STANDARD_HPP_
#define ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_STANDARD_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class TimeType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  bool using_default_policies
  >
class StepperRt
{

public:
  // required
  using time_type = TimeType;
  using state_type  = StateType;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

private:
  ::pressio::ode::StepScheme name_;
  ::pressio::ode::stepper_order_type order_;

  TimeType t_np1_  = {};
  TimeType dt_ = {};
  int32_t step_number_  = {};

  // state object to ensure the strong guarantee for handling excpetions
  state_type recovery_state_;

  // stencilStates contains:
  // bdf1: y_n
  // bdf2: y_n, y_n-1
  // cn  : y_n
  ImplicitStencilStatesContainerDyn<StateType> stencil_states_;

  // policies
  ::pressio::utils::InstanceOrReferenceWrapper<ResidualPolicyType> res_policy_;
  ::pressio::utils::InstanceOrReferenceWrapper<JacobianPolicyType> jac_policy_;

  // stencilVelocities contains:
  // bdf1,2: nothing
  // cn:     f(y_n,t_n) and f(y_np1, t_np1)
  mutable ImplicitStencilVelocitiesContainerDyn<ResidualType> stencil_velocities_;

public:
  StepperRt() = delete;
  StepperRt(const StepperRt & other)  = default;
  StepperRt & operator=(const StepperRt & other) = delete;
  StepperRt(StepperRt && other)  = default;
  StepperRt & operator=(StepperRt && other) = delete;
  ~StepperRt() = default;


  // *** BDF1 ***//
  template<
    class SystemType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::BDF1, const SystemType & systemObj)
    : name_(StepScheme::BDF1),
      order_(1),
      recovery_state_{systemObj.createState()},
      stencil_states_{systemObj.createState()},
      res_policy_{ResidualPolicyType{systemObj}},
      jac_policy_{JacobianPolicyType{systemObj}}
  {}

  template<
    class _ResidualPolicyType = ResidualPolicyType,
    class _JacobianPolicyType = JacobianPolicyType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<!_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::BDF1,
	    _ResidualPolicyType && resPolicyObj,
	    _JacobianPolicyType && jacPolicyObj)
    : name_(StepScheme::BDF1),
      order_(1),
      recovery_state_{resPolicyObj.createState()},
      stencil_states_{resPolicyObj.createState()},
      res_policy_{resPolicyObj},
      jac_policy_{jacPolicyObj}
  {}

  // *** BDF2 ***//
  template<
    class SystemType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::BDF2, const SystemType & systemObj)
    : name_(StepScheme::BDF2),
      order_(2),
      recovery_state_{systemObj.createState()},
      stencil_states_{systemObj.createState(),
		      systemObj.createState()},
      res_policy_{ResidualPolicyType{systemObj}},
      jac_policy_{JacobianPolicyType{systemObj}}
  {}

  template<
    class _ResidualPolicyType = ResidualPolicyType,
    class _JacobianPolicyType = JacobianPolicyType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<!_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::BDF2,
	    _ResidualPolicyType && resPolicyObj,
	    _JacobianPolicyType && jacPolicyObj)
    : name_(StepScheme::BDF2),
      order_(2),
      recovery_state_{resPolicyObj.createState()},
      stencil_states_{resPolicyObj.createState(),
		      resPolicyObj.createState()},
      res_policy_{resPolicyObj},
      jac_policy_{jacPolicyObj}
  {}

  // *** CN ***//
  template<
    class SystemType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::CrankNicolson, const SystemType & systemObj)
    : name_(StepScheme::CrankNicolson),
      order_(2),
      recovery_state_{systemObj.createState()},
      stencil_states_{systemObj.createState()},
      res_policy_{ResidualPolicyType{systemObj}},
      jac_policy_{JacobianPolicyType{systemObj}},
      stencil_velocities_{res_policy_.get().create(),
			  res_policy_.get().create()}
  {}

  template<
    class _ResidualPolicyType = ResidualPolicyType,
    class _JacobianPolicyType = JacobianPolicyType,
    bool _using_default_policies = using_default_policies,
    ::pressio::mpl::enable_if_t<!_using_default_policies, int > = 0
    >
  StepperRt(::pressio::ode::CrankNicolson,
	    _ResidualPolicyType && resPolicyObj,
	    _JacobianPolicyType && jacPolicyObj)
    : name_(StepScheme::CrankNicolson),
      order_(2),
      recovery_state_{resPolicyObj.createState()},
      stencil_states_{resPolicyObj.createState()},
      res_policy_{resPolicyObj},
      jac_policy_{jacPolicyObj},
      stencil_velocities_{res_policy_.get().create(),
			  res_policy_.get().create()}
  {}

public:
  ::pressio::ode::stepper_order_type order() const{
    return order_;
  }

  template<typename solver_type, typename ...Args>
  void operator()(state_type & odeState,
		  const TimeType &  currentTime,
		  const TimeType &  dt,
		  const int32_t & stepNumber,
		  solver_type & solver,
		  Args&& ...args)
  {
    PRESSIOLOG_DEBUG("implicit stepper: do step");
    auto dummyGuesser =
      [](const int32_t &, const TimeType &, state_type &){ /*no op*/ };

    if (name_==::pressio::ode::StepScheme::BDF1){
      doStepImpl(::pressio::ode::BDF1(),
		 odeState, currentTime, dt, stepNumber,
		 dummyGuesser, solver,
		 std::forward<Args>(args)...);
    }
    else if (name_==::pressio::ode::StepScheme::BDF2){
      doStepImpl(::pressio::ode::BDF2(),
		 odeState, currentTime, dt, stepNumber,
		 dummyGuesser, solver,
		 std::forward<Args>(args)...);
    }
    else if (name_==::pressio::ode::StepScheme::CrankNicolson){
      doStepImpl(::pressio::ode::CrankNicolson(),
		 odeState, currentTime, dt, stepNumber,
		 dummyGuesser, solver,
		 std::forward<Args>(args)...);
    }
  }

  // overload for when we have a guesser callback
  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void operator()(state_type & odeState,
		  const TimeType &  currentTime,
		  const TimeType &  dt,
		  const int32_t & stepNumber,
		  guess_callback_t && guesserCb,
		  solver_type & solver,
		  Args&& ...args)
  {
    PRESSIOLOG_DEBUG("implicit stepper: do step with callback to state guesser");

    if (name_==::pressio::ode::StepScheme::BDF1){
      doStepImpl(::pressio::ode::BDF1(),
		 odeState, currentTime, dt, stepNumber,
		 std::forward<guess_callback_t>(guesserCb), solver,
		 std::forward<Args>(args)...);
    }
    else if (name_==::pressio::ode::StepScheme::BDF2){
      doStepImpl(::pressio::ode::BDF2(),
		 odeState, currentTime, dt, stepNumber,
		 std::forward<guess_callback_t>(guesserCb), solver,
		 std::forward<Args>(args)...);
    }
    else if (name_==::pressio::ode::StepScheme::CrankNicolson){
      doStepImpl(::pressio::ode::CrankNicolson(),
		 odeState, currentTime, dt, stepNumber,
		 std::forward<guess_callback_t>(guesserCb), solver,
		 std::forward<Args>(args)...);
    }
  }

  auto createState() const{
    return res_policy_.get().createState();
  }

  ResidualType createResidual() const{
    return res_policy_.get().create();
  }

  JacobianType createJacobian() const{
    return jac_policy_.get().create();
  }

  void residual(const StateType & odeState, ResidualType & R) const
  {
    res_policy_.get()(name_, odeState, stencil_states_,
		      stencil_velocities_, t_np1_, dt_,
		      step_number_, R);
  }

  void jacobian(const StateType & odeState, JacobianType & J) const
  {
    jac_policy_.get()(name_, odeState, stencil_states_,
		      t_np1_, dt_, step_number_, J);
  }

private:
  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void doStepImpl(::pressio::ode::BDF1,
		  state_type & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const int32_t & stepNumber,
		  guess_callback_t && guesserCb,
		  solver_type & solver,
		  Args&& ...args)
  {

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeState.

      for bdf1, this is done by solving:  y_n+1 = y_n + dt*f(t_n+1, y_n+1)
     */

    dt_ = dt;
    t_np1_ = currentTime + dt_;
    step_number_ = stepNumber;

    // copy current solution into y_n
    auto & odeState_n = stencil_states_(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeState);

    // if provided, callback to provide a guess for the odeState
    guesserCb(stepNumber, t_np1_, odeState);

    try{
      solver.solve(*this, odeState, std::forward<Args>(args)...);
    }
    catch (::pressio::eh::NonlinearSolveFailure const & e)
    {
      // the state before attempting solution was stored in y_n-1,
      // so revert odeState to that
      auto & rollBackState = stencil_states_(ode::n());
      ::pressio::ops::deep_copy(odeState, rollBackState);

      // now throw
      throw ::pressio::eh::TimeStepFailure();
    }
  }

  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void doStepImpl(::pressio::ode::BDF2,
		  state_type & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const int32_t & stepNumber,
		  guess_callback_t && guesserCb,
		  solver_type & solver,
		  Args&& ...args)
  {

    // static_assert(::pressio::ode::legitimate_solver_for_implicit_stepper<
    //   solver_type, decltype(*this), state_type>::value,
    //   "Invalid solver for BDF2 stepper");

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeState.

      for bdf1, this is done by solving:  y_n+1 = y_n + dt*f(t_n+1, y_n+1)
     */

    dt_ = dt;
    t_np1_ = currentTime + dt;
    step_number_ = stepNumber;

    // first step, use auxiliary stepper
    if (stepNumber == 1){
        using aux_type = StepperRt<
	  TimeType, StateType, ResidualType, JacobianType,
	  const ResidualPolicyType &, const JacobianPolicyType &, false>;

	aux_type auxStepper(::pressio::ode::BDF1(),
			    res_policy_.get(), jac_policy_.get());

	// step ==1 means that we are going from y_0 to y_1
	// stencil_states_(0) now holds y_0
	::pressio::ops::deep_copy(stencil_states_(ode::n()), odeState);
	::pressio::ops::deep_copy(stencil_states_(ode::nMinusOne()), odeState);
	auxStepper(odeState, currentTime, dt,
		   stepNumber, solver,
		   std::forward<Args>(args)...);
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

      auto & odeState_n   = stencil_states_(ode::n());
      auto & odeState_nm1 = stencil_states_(ode::nMinusOne());
      ::pressio::ops::deep_copy(recovery_state_, odeState_nm1);
      ::pressio::ops::deep_copy(odeState_nm1, odeState_n);
      ::pressio::ops::deep_copy(odeState_n, odeState);

      try{
	guesserCb(stepNumber, t_np1_, odeState);
	solver.solve(*this, odeState, std::forward<Args>(args)...);
      }
      catch (::pressio::eh::NonlinearSolveFailure const & e)
	{
	  ::pressio::ops::deep_copy(odeState, odeState_n);
	  ::pressio::ops::deep_copy(odeState_n, odeState_nm1);
	  ::pressio::ops::deep_copy(odeState_nm1, recovery_state_);

	  // now throw
	  throw ::pressio::eh::TimeStepFailure();
	}
    }
  }

  template<typename solver_type, typename guess_callback_t, typename ...Args>
  void doStepImpl(::pressio::ode::CrankNicolson,
		  state_type & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const int32_t & stepNumber,
		  guess_callback_t && guesserCb,
		  solver_type & solver,
		  Args&& ...args)
  {

    /*
      upon entering this, we are at time step = stepNumber.
      The current step starts at time = currentTime and
      need to use timestep size = dt.
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeState.

      for CN, this is done by solving:
      y_n+1 = y_n + 0.5*dt*[ f(t_n+1, y_n+1) + f(t_n, y_n) ]
    */

    dt_ = dt;
    t_np1_ = currentTime + dt_;
    step_number_ = stepNumber;

    // current solution becomes y_n
    auto & odeState_n = stencil_states_(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeState);

    // if provided, callback to provide a guess for the odeState
    guesserCb(stepNumber, t_np1_, odeState);

    try{
      solver.solve(*this, odeState, std::forward<Args>(args)...);
    }
    catch (::pressio::eh::NonlinearSolveFailure const & e)
    {
      auto & rollBackState = stencil_states_(ode::n());
      ::pressio::ops::deep_copy(odeState, rollBackState);
      throw ::pressio::eh::TimeStepFailure();
    }
  }
};

}}} // end namespace pressio::ode::implicitmethods
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_STANDARD_HPP_
