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
  typename jacobian_policy_t
  >
class StepperBDF2
{
public:
  // these need to be here because are detected by solver
  using scalar_type = scalar_t;
  using state_type  = state_t;
  using residual_type = residual_t;
  using jacobian_type = jacobian_t;

  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  static constexpr types::stepper_order_t order_value = 2;
  static constexpr std::size_t numAuxStates = 2;
  using tag_name = ::pressio::ode::implicitmethods::BDF2;

  using aux_states_t =
    ::pressio::ode::AuxStatesManager<state_type, numAuxStates>;
  using standard_res_policy_t =
    ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_type, residual_type>;
  using standard_jac_policy_t =
    ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_type, jacobian_type>;

private:
  // auxiliary stepper
  aux_stepper_t & auxStepper_;

  scalar_t t_  = {};
  scalar_t dt_ = {};
  types::step_t step_  = {};
  std::reference_wrapper<const system_type> systemObj_;
  aux_states_t auxStates_;

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

  // note that here residual_policy_t can be a reference already
  // so we don't need to specify & in argument to constructor
  StepperBDF2(const state_type & state,
	      const system_type & systemObj,
	      const mpl::remove_cvref_t<residual_policy_t> & resPolicyObj,
	      const mpl::remove_cvref_t<jacobian_policy_t> & jacPolicyObj,
	      aux_stepper_t & auxStepper)
    : auxStepper_{auxStepper},
      systemObj_{systemObj},
      auxStates_{state},
      recoveryState_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
  {}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = residual_policy_t,
    typename T2 = jacobian_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, standard_res_policy_t>::value and
      mpl::is_same<T2, standard_jac_policy_t>::value,
      int> = 0
    >
  StepperBDF2(const state_type & state,
	      const system_type & systemObj,
	      aux_stepper_t & auxStepper)
    : auxStepper_{auxStepper},
      systemObj_{systemObj},
      auxStates_{state},
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
	      const scalar_t &  t,
	      const scalar_t &  dt,
	      const types::step_t & step,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("bdf2 stepper: do step");
    static_assert(::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF2 stepper");

    auto dummyGuesser =
      [](const types::step_t &, const scalar_t &, state_type &)
      { /*no op*/ };

    this->doStep(odeState, t, dt, step, solver, dummyGuesser);
  }

  // overload for when we have a guesser callback
  template<typename solver_type, typename guess_callback_t>
  void doStep(state_type & odeState,
	      const scalar_t &  t,
	      const scalar_t &  dt,
	      const types::step_t & step,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    PRESSIOLOG_DEBUG("bdf2 stepper: do step with callback to state guesser");

    static_assert(::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), state_type>::value,
      "Invalid solver for BDF2 stepper");

    using nm1 = ode::nMinusOne;
    using nm2 = ode::nMinusTwo;

    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    // first step, use auxiliary stepper
    if (step == 1){
      // step ==1 means that we are going from y_0 to y_1
      // auxStates_(0) now holds y_0
      ::pressio::ops::deep_copy(this->auxStates_(nm1()), odeState);
      auxStepper_.doStep(odeState, t, dt, step, solver);
    }
    if (step >= 2)
    {
      // step == 2 means that we are going from y_1 to y_2, so:
      //		y_n-2 = y_0 and y_n-1 = y_1
      // step == 3 means that we are going from y_2 to y_3, so:
      //		y_n-2 = y_1 and y_n-1 = y_2

      auto & odeState_nm1 = this->auxStates_(nm1());
      auto & odeState_nm2 = this->auxStates_(nm2());
      ::pressio::ops::deep_copy(recoveryState_, odeState_nm2);
      ::pressio::ops::deep_copy(odeState_nm2, odeState_nm1);
      ::pressio::ops::deep_copy(odeState_nm1, odeState);

      try{
	guesserCb(step, t, odeState);
	solver.solve(*this, odeState);
      }
      catch (::pressio::eh::nonlinear_solve_failure const & e)
	{
	  // the state before attempting solution was stored in y_n-1,
	  // so revert odeState to that
	  ::pressio::ops::deep_copy(odeState, odeState_nm1);

	  // copy y_n-2 into y_n-1
	  ::pressio::ops::deep_copy(odeState_nm1, odeState_nm2);
	  //copy safestate to y_n-2
	  ::pressio::ops::deep_copy(odeState_nm2, recoveryState_);

	  // now throw
	  throw ::pressio::eh::time_step_failure();
	}
    }
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
      <tag_name>(odeState,  this->auxStates_,
		 this->systemObj_.get(),
		 this->t_, this->dt_, this->step_, R);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const
  {
    jacPolicy_.get().template compute<
      tag_name>(odeState, this->auxStates_,
		this->systemObj_.get(),
                this->t_, this->dt_, this->step_, J);
  }
};

}}} // end namespace pressio::ode::implicitmethods
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_BDF2_IMPL_HPP_
