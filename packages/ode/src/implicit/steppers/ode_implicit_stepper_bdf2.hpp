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

#include "ode_implicit_stepper_traits_bdf2.hpp"
#include "ode_implicit_stepper_base.hpp"

namespace pressio{ namespace ode{

template<
  typename ode_state_type,
  typename ode_residual_type,
  typename ode_jacobian_type,
  typename system_type,
  typename ... Args
  >
class ImplicitStepper<
  ImplicitEnum::BDF2,
  ode_state_type,
  ode_residual_type,
  ode_jacobian_type,
  system_type,
  Args...
  >
  : public ImplicitStepperBase<
  ImplicitStepper<
    ImplicitEnum::BDF2,
    ode_state_type,
    ode_residual_type,
    ode_jacobian_type,
    system_type, Args...>
  >
{
  using this_t	       = ImplicitStepper<ImplicitEnum::BDF2,
					 ode_state_type,
					 ode_residual_type,
					 ode_jacobian_type,
					 system_type,
					 Args...>;
  using stepper_base_t = ImplicitStepperBase<this_t>;
  friend stepper_base_t;

  using mytraits       = details::traits<this_t>;
  using standard_res_policy_t = typename mytraits::standard_res_policy_t;
  using standard_jac_policy_t = typename mytraits::standard_jac_policy_t;
  using residual_pol_t = typename mytraits::residual_policy_t;
  using jacobian_pol_t = typename mytraits::jacobian_policy_t;
  using aux_stepper_t  = typename mytraits::aux_stepper_t;
  using scalar_t       = typename mytraits::scalar_t;
  static constexpr auto my_enum = mytraits::enum_id;

  aux_stepper_t & auxStepper_;

public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  ImplicitStepper(const ode_state_type & y0,
  		  const system_type & model,
  		  const residual_pol_t & resPolicyObj,
  		  const jacobian_pol_t & jacPolicyObj,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model, resPolicyObj, jacPolicyObj},
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
  ImplicitStepper(const ode_state_type & y0,
		  const system_type & model,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model},
      auxStepper_{auxStepper}{}

  // cstr for standard jacob policies
  template <
    typename T1 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepper(const ode_state_type & y0,
  		  const system_type & model,
  		  const residual_pol_t & resPolicyObj,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model, resPolicyObj},
      auxStepper_{auxStepper}{}

public:
  template<typename solver_type>
  void operator()(ode_state_type & y,
		  const scalar_t &  t,
		  const scalar_t &  dt,
		  const types::step_t & step,
		  solver_type & solver)
  {

    auto & auxY0 = this->stateAuxStorage_.data_[0];
    auto & auxY1 = this->stateAuxStorage_.data_[1];

    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    // first step, use auxiliary stepper
    if (step == 1){
      ::pressio::containers::ops::deep_copy(y, auxY0);
      auxStepper_(y, t, dt, step, solver);
    }
    if (step == 2){
      ::pressio::containers::ops::deep_copy(y, auxY1);
      solver.solve(*this, y);
    }
    if (step >= 3){
      ::pressio::containers::ops::deep_copy(auxY1, auxY0);
      ::pressio::containers::ops::deep_copy(y, auxY1);
      solver.solve(*this, y);
    }
  }

 template<typename solver_type,typename guess_callback_t>
  void operator()(ode_state_type & y,
		  const scalar_t &  t,
		  const scalar_t &  dt,
		  const types::step_t & step,
		  solver_type & solver,
		  guess_callback_t && guesserCb){

   auto & auxY0 = this->stateAuxStorage_.data_[0];
   auto & auxY1 = this->stateAuxStorage_.data_[1];

   this->dt_ = dt;
   this->t_ = t;
   this->step_t = step;

   // first step, use auxiliary stepper
   if (step == 1){
     ::pressio::containers::ops::deep_copy(y, auxY0);
     auxStepper_(y, t, dt, step, solver);
   }
   if (step == 2){
     ::pressio::containers::ops::deep_copy(y, auxY1);
     guesserCb(step, t, y);
     solver.solve(*this, y);
   }
   if (step >= 3){
     ::pressio::containers::ops::deep_copy(auxY1, auxY0);
     ::pressio::containers::ops::deep_copy(y, auxY1);
     guesserCb(step, t, y);
     solver.solve(*this, y);
   }
 }

private:
  void residualImpl(const state_type & y, residual_type & R) const
  {
    this->residual_obj_.template operator()<
      my_enum, mytraits::numAuxStates
      >(y, R, this->stateAuxStorage_.data_,
	this->sys_.get(), this->t_, this->dt_, this->step_);
  }

  residual_type residualImpl(const state_type & y) const
  {
    return this->residual_obj_.template operator()<
      my_enum, mytraits::numAuxStates
      >(y, this->stateAuxStorage_.data_,
	this->sys_.get(), this->t_, this->dt_, this->step_);
  }

  void jacobianImpl(const state_type & y, jacobian_type & J) const
  {
    this->jacobian_obj_.template operator()<
      mytraits::enum_id
      >(y, J, this->sys_.get(), this->t_, this->dt_, this->step_);
  }

  jacobian_type jacobianImpl(const state_type & y) const
  {
    return this->jacobian_obj_.template operator()<
      mytraits::enum_id
      >(y, this->sys_.get(), this->t_, this->dt_, this->step_);
  }

};//end class

}} // end namespace pressio::ode
#endif
