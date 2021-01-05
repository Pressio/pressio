/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_arbitrary_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<
  typename scalar_t,
  typename ode_state_type,
  typename ode_residual_type,
  typename ode_jacobian_type,
  typename system_type,
  typename order_setter_t,
  typename tot_n_setter_t,
  typename residual_policy_t,
  typename jacobian_policy_t,
  bool policies_are_standard
  >
class StepperArbitrary
{
public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;

  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  // numAuxStates is the number of auxiliary states needed, so all other beside y_n
  static constexpr std::size_t numAuxStates = tot_n_setter_t::value - 1;
  using tag_name = ::pressio::ode::implicitmethods::Arbitrary;
  using stencil_states_t = StencilStatesManager<ode_state_type, numAuxStates>;

private:
  scalar_t rhsEvaluationTime_  = {};
  scalar_t dt_ = {};
  types::step_t stepNumber_  = {};
  std::reference_wrapper<const system_type> systemObj_;
  stencil_states_t stencilStates_;

  // state object to ensure the strong guarantee for handling excpetions
  ode_state_type recoveryState_;
  // policies
  ::pressio::utils::instance_or_reference_wrapper<residual_policy_t> resPolicy_;
  ::pressio::utils::instance_or_reference_wrapper<jacobian_policy_t> jacPolicy_;

public:
  StepperArbitrary() = delete;
  StepperArbitrary(const StepperArbitrary & other)  = default;
  StepperArbitrary & operator=(const StepperArbitrary & other) = delete;
  StepperArbitrary(StepperArbitrary && other)  = default;
  StepperArbitrary & operator=(StepperArbitrary && other) = delete;
  ~StepperArbitrary() = default;

  // note that here residual_policy_t can be a reference already
  // so we don't need to specify & in argument to constructor
  StepperArbitrary(const ode_state_type & state,
		   const system_type & systemObj,
		   residual_policy_t resPolicyObj,
		   jacobian_policy_t jacPolicyObj)
    : systemObj_{systemObj},
      stencilStates_{state},
      recoveryState_{state},
      resPolicy_{resPolicyObj},
      jacPolicy_{jacPolicyObj}
    {}

  // cstr for standard residual and jacob policies
  template <
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int > = 0
    >
  StepperArbitrary(const ode_state_type & state,
                   const system_type & systemObj)
    : systemObj_{systemObj},
      stencilStates_{state},
      recoveryState_{state},
      resPolicy_{},
      jacPolicy_{}
  {}

public:
  ::pressio::ode::types::stepper_order_t order() const
  {
    return order_setter_t::value;
  }

  residual_type createResidual() const
  {
    return resPolicy_.get().create(systemObj_.get());
  }

  jacobian_type createJacobian() const
  {
    return jacPolicy_.get().create(systemObj_.get());
  }

  void residual(const state_type & odeState, residual_type & R) const
  {
    resPolicy_.get().template compute
      <tag_name>(odeState, stencilStates_,
		 systemObj_.get(),
		 rhsEvaluationTime_, dt_, stepNumber_, R);
  }

  void jacobian(const state_type & odeState, jacobian_type & J) const
  {
    jacPolicy_.get().template compute
      <tag_name>(odeState, stencilStates_,
		 systemObj_.get(),
		 rhsEvaluationTime_, dt_, stepNumber_, J);
  }

  template<typename solver_type>
  void doStep(ode_state_type & odeState,
	      const scalar_t & currentTime,
	      const scalar_t & dt,
	      const types::step_t & step,
	      solver_type & solver)
  {
    PRESSIOLOG_DEBUG("arbitrary stepper: do step");

    static_assert
      (::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), ode_state_type>::value,
      "Invalid solver for Arbitrary stepper");

    dt_ = dt;
    rhsEvaluationTime_ = currentTime + dt;
    stepNumber_ = step;

    updateAuxiliaryStorage<numAuxStates>(odeState);

    try{
      solver.solve(*this, odeState);
    }
    catch (::pressio::eh::nonlinear_solve_failure const & e)
    {
      rollBackStates<numAuxStates>(odeState);
      throw ::pressio::eh::time_step_failure();
    }
  }

private:
  // one aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==1>
  updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    ::pressio::ops::deep_copy(recoveryState_, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==1>
  rollBackStates(ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, recoveryState_);
  }

  // two aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==2>
  updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    auto & y_nm1 = stencilStates_.stateAt(ode::nMinusOne());
    ::pressio::ops::deep_copy(recoveryState_, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==2>
  rollBackStates(ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    auto & y_nm1 = stencilStates_.stateAt(ode::nMinusOne());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, recoveryState_);
  }

  // three aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==3>
  updateAuxiliaryStorage(const ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    auto & y_nm1 = stencilStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = stencilStates_.stateAt(ode::nMinusTwo());
    ::pressio::ops::deep_copy(recoveryState_, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==3>
  rollBackStates(ode_state_type & odeState)
  {
    auto & y_n = stencilStates_.stateAt(ode::n());
    auto & y_nm1 = stencilStates_.stateAt(ode::nMinusOne());
    auto & y_nm2 = stencilStates_.stateAt(ode::nMinusTwo());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, recoveryState_);
  }
};

}}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
