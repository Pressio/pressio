/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_arbitrary.hpp
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

#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_ARBITRARY_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_ARBITRARY_HPP_

#include "ode_implicit_stepper_traits_arbitrary.hpp"
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
  ImplicitEnum::Arbitrary,
  ode_state_type,
  ode_residual_type,
  ode_jacobian_type,
  system_type,
  Args...
  >
  : public ImplicitStepperBase<
  ImplicitStepper<
    ImplicitEnum::Arbitrary, ode_state_type,
    ode_residual_type, ode_jacobian_type,
    system_type, Args...>
  >
{

  using this_t	       = ImplicitStepper<ImplicitEnum::Arbitrary,
					 ode_state_type,
					 ode_residual_type,
					 ode_jacobian_type,
					 system_type,
					 Args...>;
  using stepper_base_t = ImplicitStepperBase<this_t>;
  friend stepper_base_t;

  using mytraits		= details::traits<this_t>;
  using standard_res_policy_t	= typename mytraits::standard_res_policy_t;
  using standard_jac_policy_t	= typename mytraits::standard_jac_policy_t;
  using residual_pol_t		= typename mytraits::residual_policy_t;
  using jacobian_pol_t		= typename mytraits::jacobian_policy_t;
  using scalar_t		= typename mytraits::scalar_t;
  static constexpr auto my_enum = mytraits::enum_id;

public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  ImplicitStepper(const ode_state_type & stateIn0,
  		  const system_type & model,
  		  const residual_pol_t & resPolicyObj,
  		  const jacobian_pol_t & jacPolicyObj)
    : stepper_base_t{stateIn0, model, resPolicyObj, jacPolicyObj}{}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepper(const ode_state_type & stateIn0,
		  const system_type & model)
    : stepper_base_t{stateIn0, model}{}

public:
  template<typename solver_type>
  void operator()(ode_state_type & odeState,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const types::step_t & step,
		  solver_type & solver)
  {
    this->dt_ = dt;
    this->t_ = t;
    this->step_ = step;

    constexpr auto nAux = decltype(this->auxStates_)::size();
    this->updateAuxiliaryStorage<nAux>(odeState);
    solver.solve(*this, odeState);
  }

private:
  // methods to do updating on the storage of previous states
  template<std::size_t nAux, mpl::enable_if_t<nAux==1> * = nullptr>
  void updateAuxiliaryStorage(const ode_state_type & odeState){
    // copy y_n into y_n-1
    auto & auxY0 = this->auxStates_[0];
    ::pressio::containers::ops::deep_copy(odeState, auxY0);
  }

  // when we have two aux states,
  template<std::size_t nAux, mpl::enable_if_t<nAux==2> * = nullptr>
  void updateAuxiliaryStorage(const ode_state_type & odeState){
    // copy y_n-1 into y_n-2
    auto & aux1 = this->auxStates_[0];
    auto & aux2 = this->auxStates_[1];
    ::pressio::containers::ops::deep_copy(aux1, aux2);
    // copy y_n into y_n-1
    ::pressio::containers::ops::deep_copy(odeState, aux1);
  }

  // when we have three aux states,
  template<std::size_t nAux, mpl::enable_if_t<nAux==3> * = nullptr>
  void updateAuxiliaryStorage(const ode_state_type & odeState){
    auto & aux1 = this->auxStates_[0];
    auto & aux2 = this->auxStates_[1];
    auto & aux3 = this->auxStates_[2];
    ::pressio::containers::ops::deep_copy(aux2, aux3);
    ::pressio::containers::ops::deep_copy(aux1, aux2);
    ::pressio::containers::ops::deep_copy(odeState, aux1);
  }

  // when we have four aux states,
  template<std::size_t nAux, mpl::enable_if_t<nAux==4> * = nullptr>
  void updateAuxiliaryStorage(const ode_state_type & odeState){
    auto & aux1 = this->auxStates_[0];
    auto & aux2 = this->auxStates_[1];
    auto & aux3 = this->auxStates_[2];
    auto & aux4 = this->auxStates_[3];
    ::pressio::containers::ops::deep_copy(aux3, aux4);
    ::pressio::containers::ops::deep_copy(aux2, aux3);
    ::pressio::containers::ops::deep_copy(aux1, aux2);
    ::pressio::containers::ops::deep_copy(odeState, aux1);
  }

  // // when we have five aux states,
  // template<std::size_t nAux, mpl::enable_if_t<nAux==5> * = nullptr>
  // void updateAuxiliaryStorage(const ode_state_type & odeState){
  //   for (auto i=nAux-2; i>=0; --i){
  //     auto & source	  = this->auxStates_[i];
  //     auto & destination = this->auxStates_[i+1];
  //     ::pressio::containers::ops::deep_copy(source, destination);
  //   }
  //   auto & aux1 = this->auxStates_[0];
  //   ::pressio::containers::ops::deep_copy(odeState, aux1);
  // }


private:
  void residualImpl(const state_type & odeState, residual_type & R) const
  {
    this->residual_obj_.template operator()<
      mytraits::numAuxStates
      >(odeState, this->auxStates_,
       this->sys_.get(), this->t_, this->dt_, this->step_, R);
  }

  residual_type residualImpl(const state_type & odeState) const
  {
    return this->residual_obj_.operator()(odeState, this->sys_.get(),
					  this->t_, this->step_);

    // return this->residual_obj_.template operator()<
    //   mytraits::numAuxStates
    //   >(odeState, this->auxStates_,
    //    this->sys_.get(), this->t_, this->dt_, this->step_);
  }

  void jacobianImpl(const state_type & odeState, jacobian_type & J) const
  {
    this->jacobian_obj_.template operator()<
      mytraits::numAuxStates
      >(odeState, this->auxStates_,
       this->sys_.get(), this->t_, this->dt_, this->step_, J);
  }

  jacobian_type jacobianImpl(const state_type & odeState) const
  {
    return this->jacobian_obj_.operator()(odeState, this->sys_.get(),
					  this->dt_, this->step_);

    // return this->jacobian_obj_.template operator()<
    //   mytraits::numAuxStates
    //   >(odeState, this->auxStates_,
    // 	this->sys_.get(), this->t_, this->dt_, this->step_);
  }

};//end class

}} // end namespace pressio::ode
#endif
