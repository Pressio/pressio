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

#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_BDF2_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_BDF2_HPP_

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

public:
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using tag_name = ::pressio::ode::implicitmethods::BDF2;

  static constexpr types::stepper_order_t order_value = 2;
  static constexpr std::size_t numAuxStates = 2;
  using system_wrapper_t  = ::pressio::ode::impl::OdeSystemWrapper<system_type>;
  using aux_states_t = ::pressio::ode::AuxStatesManager<state_type, numAuxStates>;

  using standard_res_policy_t = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_type, system_type, residual_type>;
  using standard_jac_policy_t = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_type, system_type, jacobian_type>;

public:
  StepperBDF2() = delete;
  ~StepperBDF2() = default;
  StepperBDF2(const StepperBDF2 & other)  = delete;
  StepperBDF2 & operator=(const StepperBDF2 & other)  = delete;
  StepperBDF2(StepperBDF2 && other)  = delete;
  StepperBDF2 & operator=(StepperBDF2 && other)  = delete;

  StepperBDF2(const state_type & state,
	  const system_type & systemObj,
	  const residual_policy_t & resPolicyObj,
	  const jacobian_policy_t & jacPolicyObj,
	  aux_stepper_t & auxStepper)
    : auxStepper_{auxStepper}, 
      sys_{systemObj}, 
      auxStates_{state}, 
      residual_obj_{resPolicyObj},
      jacobian_obj_{jacPolicyObj}
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
      sys_{systemObj}, 
      auxStates_{state}
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
    static_assert(::pressio::ode::concepts::legitimate_solver_for_implicit_stepper<
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
      ::pressio::ops::deep_copy(this->auxStates_.get(nm1()), odeState);
      // advnace the ode state
      auxStepper_.doStep(odeState, t, dt, step, solver);
    }
    if (step >= 2)
    {
      // step == 2 means that we are going from y_1 to y_2, so:
      //		y_n-2 = y_0 and y_n-1 = y_1
      // step == 3 means that we are going from y_2 to y_3, so:
      //		y_n-2 = y_1 and y_n-1 = y_2

      auto & odeState_nm1 = this->auxStates_.get(nm1());
      auto & odeState_nm2 = this->auxStates_.get(nm2());
      ::pressio::ops::deep_copy(odeState_nm2, odeState_nm1);
      ::pressio::ops::deep_copy(odeState_nm1, odeState);
      solver.solve(*this, odeState);
    }
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
    static_assert(::pressio::ode::concepts::legitimate_solver_for_implicit_stepper<
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
      ::pressio::ops::deep_copy(this->auxStates_.get(nm1()), odeState);
      // advnace the ode state
      auxStepper_(odeState, t, dt, step, solver);
    }
    if (step >= 2)
    {
      // step == 2 means that we are going from y_1 to y_2, so:
      //		y_n-2 = y_0 and y_n-1 = y_1
      // step == 3 means that we are going from y_2 to y_3, so:
      //		y_n-2 = y_1 and y_n-1 = y_2

      auto & odeState_nm1 = this->auxStates_.get(nm1());
      auto & odeState_nm2 = this->auxStates_.get(nm2());
      ::pressio::ops::deep_copy(odeState_nm2, odeState_nm1);
      ::pressio::ops::deep_copy(odeState_nm1, odeState);
      guesserCb(step, t, odeState);
      solver.solve(*this, odeState);
    }
  }

  residual_t createResidual() const
  {
    return this->residual_obj_.create(sys_.get());
  }

  jacobian_t createJacobian() const
  {
    return this->jacobian_obj_.create(sys_.get());
  }

  void residual(const state_t & odeState, residual_t & R,
    ::pressio::Norm normKind, scalar_t & normValue) const
  {
    this->residual_obj_.template compute<tag_name>(
      odeState,  this->auxStates_, this->sys_.get(),
      this->t_, this->dt_, this->step_, R,
      normKind, normValue);
  }

  void jacobian(const state_t & odeState, jacobian_t & J) const
  {
    this->jacobian_obj_.template compute<
      tag_name>(odeState, this->auxStates_, this->sys_.get(), 
                this->t_, this->dt_, this->step_, J);
  }

private:
  aux_stepper_t & auxStepper_;

  scalar_t t_  = {};
  scalar_t dt_ = {};
  types::step_t step_  = {};
  system_wrapper_t sys_;
  aux_states_t auxStates_;

  // conditionally set the type of the object knowing how to compute residual
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_res_policy_t, residual_policy_t>::value,
    const residual_policy_t,
    const residual_policy_t &
    >::type residual_obj_;

  // conditionally set the type of the object knowing how to compute jacobian
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_jac_policy_t, jacobian_policy_t>::value,
    const jacobian_policy_t,
    const jacobian_policy_t &
    >::type jacobian_obj_;
};//end class

}}} // end namespace pressio::ode::implicitmethods
#endif
