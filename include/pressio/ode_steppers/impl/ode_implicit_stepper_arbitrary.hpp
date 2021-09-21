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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_IMPL_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  int n_states,
  class ScalarType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class SystemType
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
  using tag_name = ::pressio::ode::ImplicitArbitrary;
  using stencil_states_t = ImplicitStencilStatesContainerStatic<StateType, numAuxStates>;

private:
  ScalarType rhsEvaluationTime_  = {};
  ScalarType dt_ = {};
  int32_t stepNumber_  = {};

  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  stencil_states_t stencilStates_;
  // state object to ensure the strong guarantee for handling excpetions
  StateType recoveryState_;

public:
  StepperArbitrary() = delete;
  StepperArbitrary(const StepperArbitrary & other)  = default;
  StepperArbitrary & operator=(const StepperArbitrary & other) = delete;
  StepperArbitrary(StepperArbitrary && other)  = default;
  StepperArbitrary & operator=(StepperArbitrary && other) = delete;
  ~StepperArbitrary() = default;

  StepperArbitrary(const StateType & state,
		               SystemType && systemObj)
    : systemObj_(std::forward<SystemType>(systemObj)),
      stencilStates_(state), //stencilstates handles right semantics
      recoveryState_{::pressio::ops::clone(state)}
    {}

public:
  residual_type createResidual() const{
    return systemObj_.get().createDiscreteTimeResidual();
  }

  jacobian_type createJacobian() const{
    return systemObj_.get().createDiscreteTimeJacobian();
  }

  // 1 aux states
  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==1 >
  residual(const state_type & odeState, residual_type & R) const
  {
    const auto & yn = stencilStates_(ode::n());

    try{
    systemObj_.get().template discreteTimeResidual(stepNumber_, rhsEvaluationTime_,
      dt_, R, odeState, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==1 >
  jacobian(const state_type & odeState, jacobian_type & J) const
  {
    const auto & yn = stencilStates_(ode::n());
    systemObj_.get().template discreteTimeJacobian(stepNumber_, rhsEvaluationTime_,
      dt_, J, odeState, yn);
  }

  // 2 aux states
  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==2 >
  residual(const state_type & odeState, residual_type & R) const
  {
    const auto & yn = stencilStates_(ode::n());
    const auto & ynm1 = stencilStates_(ode::nMinusOne());

    try{
    systemObj_.get().template discreteTimeResidual(stepNumber_, rhsEvaluationTime_,
      dt_, R, odeState, yn, ynm1);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==2 >
  jacobian(const state_type & odeState, jacobian_type & J) const
  {
    const auto & yn = stencilStates_(ode::n());
    const auto & ynm1 = stencilStates_(ode::nMinusOne());
    systemObj_.get().template discreteTimeJacobian(stepNumber_, rhsEvaluationTime_,
      dt_, J, odeState, yn, ynm1);
  }

  // 3 aux states
  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==3 >
  residual(const state_type & odeState, residual_type & R) const
  {
    const auto & yn = stencilStates_(ode::n());
    const auto & ynm1 = stencilStates_(ode::nMinusOne());
    const auto & ynm2 = stencilStates_(ode::nMinusTwo());

    try{
    systemObj_.get().template discreteTimeResidual(stepNumber_, rhsEvaluationTime_,
      dt_, R, odeState, yn, ynm1, ynm2);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template< std::size_t _numAuxStates = numAuxStates>
  mpl::enable_if_t< _numAuxStates==3 >
  jacobian(const state_type & odeState, jacobian_type & J) const
  {
    const auto & yn = stencilStates_(ode::n());
    const auto & ynm1 = stencilStates_(ode::nMinusOne());
    const auto & ynm2 = stencilStates_(ode::nMinusTwo());
    systemObj_.get().template discreteTimeJacobian(stepNumber_, rhsEvaluationTime_,
      dt_, J, odeState, yn, ynm1, ynm2);
  }


  template<typename SolverType, class StepCountType, typename ...Args>
  void operator()(StateType & odeState,
	      const ScalarType & currentTime,
	      const ScalarType & dt,
	      const StepCountType & step,
	      SolverType & solver,
	      Args&& ...args)

  {
    PRESSIOLOG_DEBUG("arbitrary stepper: do step");
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
