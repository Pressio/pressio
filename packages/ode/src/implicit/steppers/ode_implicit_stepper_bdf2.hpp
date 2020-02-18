/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_bdf2.hpp
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
  typename ode_state_type,
  typename ode_residual_type,
  typename ode_jacobian_type,
  typename system_type,
  typename ... Args
  >
class Stepper<
  ::pressio::ode::implicitmethods::BDF2,
  ode_state_type,
  ode_residual_type,
  ode_jacobian_type,
  system_type,
  Args...
  >
  : public StepperBase<
  Stepper<
    ::pressio::ode::implicitmethods::BDF2,
    ode_state_type,
    ode_residual_type,
    ode_jacobian_type,
    system_type, Args...>
  >
{
  using this_t = Stepper<::pressio::ode::implicitmethods::BDF2,
			  ode_state_type,
			  ode_residual_type,
			  ode_jacobian_type,
			  system_type,
			  Args...>;
  using stepper_base_t = StepperBase<this_t>;
  friend stepper_base_t;
  using typename stepper_base_t::aux_states_t;

  using mytraits       = details::traits<this_t>;
  using standard_res_policy_t = typename mytraits::standard_res_policy_t;
  using standard_jac_policy_t = typename mytraits::standard_jac_policy_t;
  using residual_pol_t = typename mytraits::residual_policy_t;
  using jacobian_pol_t = typename mytraits::jacobian_policy_t;
  using aux_stepper_t  = typename mytraits::aux_stepper_t;
  using scalar_t       = typename mytraits::scalar_t;
  using tag_name       = typename mytraits::tag_name;

  aux_stepper_t & auxStepper_;

public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;

public:
  Stepper() = delete;
  ~Stepper() = default;

  // copy cnstr
  Stepper(const Stepper & other)  = delete;
  // copy assignment
  Stepper & operator=(const Stepper & other)  = delete;
  // move cnstr
  Stepper(Stepper && other)  = delete;
  // move assign
  Stepper & operator=(Stepper && other)  = delete;

  Stepper(const ode_state_type & stateIn0,
	  const system_type & model,
	  const residual_pol_t & resPolicyObj,
	  const jacobian_pol_t & jacPolicyObj,
	  aux_stepper_t & auxStepper)
    : stepper_base_t{stateIn0, model, resPolicyObj, jacPolicyObj},
      auxStepper_{auxStepper}{}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  Stepper(const ode_state_type & stateIn0,
	  const system_type & model,
	  aux_stepper_t & auxStepper)
    : stepper_base_t{stateIn0, model},
      auxStepper_{auxStepper}{}

  // cstr for standard jacob policies
  template <
    typename T1 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, jacobian_pol_t>::value
      > * = nullptr
    >
  Stepper(const ode_state_type & stateIn0,
	  const system_type & model,
	  const residual_pol_t & resPolicyObj,
	  aux_stepper_t & auxStepper)
    : stepper_base_t{stateIn0, model, resPolicyObj},
      auxStepper_{auxStepper}{}

private:
  // enable when auxiliary stepper is implicit too
  template<typename solver_type>
  void doStep(ode_state_type & odeState,
	      const scalar_t &  t,
	      const scalar_t &  dt,
	      const types::step_t & step,
	      solver_type & solver)
  {
    using nm1 = ode::nMinusOne;
    using nm2 = ode::nMinusTwo;

    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    // first step, use auxiliary stepper
    if (step == 1){
      // step ==1 means that we are going from y_0 to y_1
      // auxStates_(0) now holds y_0
      ::pressio::containers::ops::deep_copy(odeState, this->auxStates_.get(nm1()));
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
      ::pressio::containers::ops::deep_copy(odeState_nm1, odeState_nm2);
      ::pressio::containers::ops::deep_copy(odeState, odeState_nm1);
      solver.solve(*this, odeState);
    }
  }


  // enable when auxiliary stepper is implicit too
  // overload for when we have a guesser callback
  template<typename solver_type, typename guess_callback_t>
  void doStep(ode_state_type & odeState,
	      const scalar_t &  t,
	      const scalar_t &  dt,
	      const types::step_t & step,
	      solver_type & solver,
	      guess_callback_t && guesserCb)
  {
    using nm1 = ode::nMinusOne;
    using nm2 = ode::nMinusTwo;

    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    // first step, use auxiliary stepper
    if (step == 1){
      // step ==1 means that we are going from y_0 to y_1
      // auxStates_(0) now holds y_0
      ::pressio::containers::ops::deep_copy(odeState, this->auxStates_.get(nm1()));
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
      ::pressio::containers::ops::deep_copy(odeState_nm1, odeState_nm2);
      ::pressio::containers::ops::deep_copy(odeState, odeState_nm1);
      guesserCb(step, t, odeState);
      solver.solve(*this, odeState);
    }
  }

  void residualImpl(const state_type & odeState, residual_type & R) const
  {
    this->residual_obj_.template operator()<
      tag_name>(odeState, this->auxStates_,
		this->sys_.get(), this->t_, this->dt_, this->step_, R);
  }

  residual_type residualImpl(const state_type & odeState) const
  {
    return this->residual_obj_.template operator()<
      tag_name>(odeState, this->auxStates_,
		this->sys_.get(), this->t_, this->dt_, this->step_);
  }

  void jacobianImpl(const state_type & odeState, jacobian_type & J) const
  {
    this->jacobian_obj_.template operator()<
      tag_name
      >(odeState, this->sys_.get(), this->t_, this->dt_, this->step_, J);
  }

  jacobian_type jacobianImpl(const state_type & odeState) const
  {
    return this->jacobian_obj_.template operator()<
      tag_name
      >(odeState, this->sys_.get(), this->t_, this->dt_, this->step_);
  }

};//end class

}}} // end namespace pressio::ode::implicitmethods
#endif
