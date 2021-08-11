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

namespace pressio{ namespace ode{ namespace impl{

template<
  int order_in,
  int n_states,
  class ScalarType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class SystemType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  bool policies_are_standard
  >
class StepperArbitrary
{
public:
  // these need to be here because are detected by solver
  using scalar_type	= ScalarType;
  using state_type	= StateType;
  using residual_type	= ResidualType;
  using jacobian_type	= JacobianType;

  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  // numAuxStates is the number of auxiliary states needed, so all other beside y_n
  static constexpr std::size_t numAuxStates = n_states - 1;
  using tag_name = ::pressio::ode::implicitmethods::Arbitrary;
  using stencil_states_t = ImplicitStencilStatesContainer<StateType, numAuxStates>;

private:
  ScalarType rhsEvaluationTime_  = {};
  ScalarType dt_ = {};
  ::pressio::ode::step_count_type stepNumber_  = {};
  std::reference_wrapper<const SystemType> systemObj_;
  stencil_states_t stencilStates_;

  // state object to ensure the strong guarantee for handling excpetions
  StateType recoveryState_;
  // policies
  ::pressio::utils::InstanceOrReferenceWrapper<ResidualPolicyType> resPolicy_;
  ::pressio::utils::InstanceOrReferenceWrapper<JacobianPolicyType> jacPolicy_;

public:
  StepperArbitrary() = delete;
  StepperArbitrary(const StepperArbitrary & other)  = default;
  StepperArbitrary & operator=(const StepperArbitrary & other) = delete;
  StepperArbitrary(StepperArbitrary && other)  = default;
  StepperArbitrary & operator=(StepperArbitrary && other) = delete;
  ~StepperArbitrary() = default;

  // note that here ResidualPolicyType can be a reference already
  // so we don't need to specify & in argument to constructor
  StepperArbitrary(const StateType & state,
		   const SystemType & systemObj,
		   ResidualPolicyType && resPolicyObj,
		   JacobianPolicyType && jacPolicyObj)
    : systemObj_{systemObj},
      stencilStates_(state), //stencilstates handles right semantics
      recoveryState_{::pressio::ops::clone(state)},
      resPolicy_{std::forward<ResidualPolicyType>(resPolicyObj)},
      jacPolicy_{std::forward<JacobianPolicyType>(jacPolicyObj)}
    {}

  // cstr for standard residual and jacob policies
  template <
    bool _policies_are_standard = policies_are_standard,
    ::pressio::mpl::enable_if_t<_policies_are_standard, int > = 0
    >
  StepperArbitrary(const StateType & state,
                   const SystemType & systemObj)
    : systemObj_{systemObj},
      stencilStates_(state), //stencilstates handles right semantics
      recoveryState_{::pressio::ops::clone(state)},
      resPolicy_{},
      jacPolicy_{}
  {}

public:
  ::pressio::ode::stepper_order_type order() const
  {
    return order_in;
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

  template<typename solver_type, typename ...Args>
  void doStep(StateType & odeState,
	      const ScalarType & currentTime,
	      const ScalarType & dt,
	      const ::pressio::ode::step_count_type & step,
	      solver_type & solver,
	      Args&& ...args)

  {
    PRESSIOLOG_DEBUG("arbitrary stepper: do step");

    static_assert
      (::pressio::ode::legitimate_solver_for_implicit_stepper<
      solver_type, decltype(*this), StateType>::value,
      "Invalid solver for Arbitrary stepper");

    dt_ = dt;
    rhsEvaluationTime_ = currentTime + dt;
    stepNumber_ = step;

    updateAuxiliaryStorage<numAuxStates>(odeState);

    try{
      solver.solve(*this, odeState, std::forward<Args>(args)...);
    }
    catch (::pressio::eh::NonlinearSolveFailure const & e)
    {
      rollBackStates<numAuxStates>(odeState);
      throw ::pressio::eh::TimeStepFailure();
    }
  }

private:
  // one aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==1>
  updateAuxiliaryStorage(const StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    ::pressio::ops::deep_copy(recoveryState_, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==1>
  rollBackStates(StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, recoveryState_);
  }

  // two aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==2>
  updateAuxiliaryStorage(const StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    auto & y_nm1 = stencilStates_(ode::nMinusOne());
    ::pressio::ops::deep_copy(recoveryState_, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==2>
  rollBackStates(StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    auto & y_nm1 = stencilStates_(ode::nMinusOne());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, recoveryState_);
  }

  // three aux states
  template<std::size_t nAux>
  mpl::enable_if_t<nAux==3>
  updateAuxiliaryStorage(const StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    auto & y_nm1 = stencilStates_(ode::nMinusOne());
    auto & y_nm2 = stencilStates_(ode::nMinusTwo());
    ::pressio::ops::deep_copy(recoveryState_, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  mpl::enable_if_t<nAux==3>
  rollBackStates(StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    auto & y_nm1 = stencilStates_(ode::nMinusOne());
    auto & y_nm2 = stencilStates_(ode::nMinusTwo());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, y_nm1);
    ::pressio::ops::deep_copy(y_nm1, y_nm2);
    ::pressio::ops::deep_copy(y_nm2, recoveryState_);
  }
};

}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
