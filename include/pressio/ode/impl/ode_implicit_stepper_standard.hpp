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
  class IndVarType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class ResidualJacobianPolicyType
  >
class ImplicitStepperStandardImpl
{

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type  = StateType;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

private:
  ::pressio::ode::StepScheme name_;

  IndVarType t_np1_  = {};
  IndVarType dt_ = {};
  int32_t step_number_  = {};

  // state object to ensure the strong guarantee for handling excpetions
  state_type recovery_state_;

  // stencilStates contains:
  // for bdf1: y_n
  // for bdf2: y_n, y_n-1
  // for cn  : y_n
  ImplicitStencilStatesDynamicContainer<StateType> stencil_states_;

  ::pressio::utils::InstanceOrReferenceWrapper<ResidualJacobianPolicyType> rj_policy_;

  // stencilRightHandSide contains:
  // for bdf1,2: nothing
  // for cn:  f(y_n,t_n) and f(y_np1, t_np1)
  mutable ImplicitStencilRightHandSideDynamicContainer<ResidualType> stencil_rhs_;

public:
  ImplicitStepperStandardImpl() = delete;
  ImplicitStepperStandardImpl(const ImplicitStepperStandardImpl & other)  = default;
  ImplicitStepperStandardImpl & operator=(const ImplicitStepperStandardImpl & other) = delete;
  ~ImplicitStepperStandardImpl() = default;

  // *** BDF1 ***//
  ImplicitStepperStandardImpl(::pressio::ode::BDF1,
			      ResidualJacobianPolicyType && rjPolicyObj)
    : name_(StepScheme::BDF1),
      recovery_state_{rjPolicyObj.createState()},
      stencil_states_{rjPolicyObj.createState()},
      rj_policy_(std::forward<ResidualJacobianPolicyType>(rjPolicyObj))
  {}

  // *** BDF2 ***//
  ImplicitStepperStandardImpl(::pressio::ode::BDF2,
			      ResidualJacobianPolicyType && rjPolicyObj)
    : name_(StepScheme::BDF2),
      recovery_state_{rjPolicyObj.createState()},
      stencil_states_{rjPolicyObj.createState(),
		      rjPolicyObj.createState()},
      rj_policy_(std::forward<ResidualJacobianPolicyType>(rjPolicyObj))
  {}

  // *** CN ***//
  ImplicitStepperStandardImpl(::pressio::ode::CrankNicolson,
			      ResidualJacobianPolicyType && rjPolicyObj)
    : name_(StepScheme::CrankNicolson),
      recovery_state_{rjPolicyObj.createState()},
      stencil_states_{rjPolicyObj.createState()},
      rj_policy_(std::forward<ResidualJacobianPolicyType>(rjPolicyObj)),
      stencil_rhs_{rj_policy_.get().createResidual(),
                   rj_policy_.get().createResidual()}
  {}

public:
  template<class SolverType, class ...Args>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount stepNumber,
		  const ::pressio::ode::StepSize<independent_variable_type> & stepSize,
		  SolverType & solver,
		  Args && ...args)
  {
    PRESSIOLOG_DEBUG("implicit stepper: do step");

    if (name_==::pressio::ode::StepScheme::BDF1){
      doStepImpl(::pressio::ode::BDF1(),
		 odeState, stepStartVal.get(), stepSize.get(),
		 stepNumber.get(), solver,
		 std::forward<Args>(args)...);
    }

    else if (name_==::pressio::ode::StepScheme::BDF2){
      doStepImpl(::pressio::ode::BDF2(),
		 odeState, stepStartVal.get(), stepSize.get(),
		 stepNumber.get(), solver,
		 std::forward<Args>(args)...);
    }

    else if (name_==::pressio::ode::StepScheme::CrankNicolson){
      doStepImpl(::pressio::ode::CrankNicolson(),
		 odeState, stepStartVal.get(), stepSize.get(),
		 stepNumber.get(), solver,
		 std::forward<Args>(args)...);
    }
  }

  auto createState() const{
    return rj_policy_.get().createState();
  }

  ResidualType createResidual() const{
    return rj_policy_.get().createResidual();
  }

  JacobianType createJacobian() const{
    return rj_policy_.get().createJacobian();
  }

  void residualAndJacobian(const StateType & odeState,
			   ResidualType & R,
			   JacobianType & J,
			   bool recomputeJacobian) const
  {

    rj_policy_.get()(name_, odeState, stencil_states_,
		     stencil_rhs_,
		     ::pressio::ode::StepEndAt<IndVarType>(t_np1_),
		     ::pressio::ode::StepCount(step_number_),
		     ::pressio::ode::StepSize<IndVarType>(dt_),
		     R, J, recomputeJacobian);
  }

private:
  template<class solver_type, class ...Args>
  void doStepImpl(::pressio::ode::BDF1,
		  state_type & odeState,
		  const IndVarType & currentTime,
		  const IndVarType & dt,
		  const int32_t & stepNumber,
		  solver_type & solver,
		  Args&& ...args)
  {

    /*
      here, we are at step = stepNumber.
      The current step starts at time = currentTime and
      we need to use timestep size = dt.

      bdf1 predicts next state y_n+1 by solving:
          R = y_n+1 - y_n - dt*f(t_n+1, y_n+1) = 0
      predict/compute the solution at the next step: stepNumber+1
      and store that into odeState.
     */

    dt_ = dt;
    t_np1_ = currentTime + dt_;
    step_number_ = stepNumber;

    // copy current solution into y_n
    auto & odeState_n = stencil_states_(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeState);

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

  template<class solver_type, class ...Args>
  void doStepImpl(::pressio::ode::BDF2,
		  state_type & odeState,
		  const IndVarType & currentTime,
		  const IndVarType & dt,
		  const int32_t & stepNumber,
		  solver_type & solver,
		  Args&& ...args)
  {

    dt_ = dt;
    t_np1_ = currentTime + dt;
    step_number_ = stepNumber;

    // first step, use auxiliary stepper
    if (stepNumber == 1){
        using aux_type = ImplicitStepperStandardImpl<
	  IndVarType, StateType, ResidualType, JacobianType,
	  const ResidualJacobianPolicyType &>;

	aux_type auxStepper(::pressio::ode::BDF1(), rj_policy_.get());

	// step ==1 means that we are going from y_0 to y_1
	// copy y_0 into stencil_states_(0)
	::pressio::ops::deep_copy(stencil_states_(ode::n()), odeState);
	::pressio::ops::deep_copy(stencil_states_(ode::nMinusOne()), odeState);
	auxStepper(odeState, StepStartAt<independent_variable_type>(currentTime),
		   StepCount(stepNumber),
		   StepSize<independent_variable_type>(dt),
		   solver, std::forward<Args>(args)...);
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

  template<class solver_type, class ...Args>
  void doStepImpl(::pressio::ode::CrankNicolson,
		  state_type & odeState,
		  const IndVarType & currentTime,
		  const IndVarType & dt,
		  const int32_t & stepNumber,
		  solver_type & solver,
		  Args&& ...args)
  {

    /*
      y_n+1 = y_n + 0.5*dt*[ f(t_n+1, y_n+1) + f(t_n, y_n) ]
    */

    dt_ = dt;
    t_np1_ = currentTime + dt_;
    step_number_ = stepNumber;

    // current solution becomes y_n
    auto & odeState_n = stencil_states_(ode::n());
    ::pressio::ops::deep_copy(odeState_n, odeState);

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
