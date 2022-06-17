/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_WITH_MASS_MATRIX_HPP_
#define ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_WITH_MASS_MATRIX_HPP_

#include <array>

namespace pressio{ namespace ode{ namespace impl{

// this class is not meant for direct instantiation.
// One needs to use the public create_* functions because
// templates are handled and passed properly there.
template<
  class MassMatrixType,
  class StateType,
  class IndVarType,
  class SystemType,
  class RightHandSideType
  >
class ExplicitStepperWithMassMatrix
{
public:
  using independent_variable_type  = IndVarType;
  using state_type  = StateType;

private:
  StepScheme name_;
  const stepper_order_type order_;
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  std::vector<RightHandSideType> rhsInstances_;
  StateType auxiliaryState_;
  MassMatrixType massMatrix_;

public:
  ExplicitStepperWithMassMatrix() = delete;
  ExplicitStepperWithMassMatrix(const ExplicitStepperWithMassMatrix &) = default;
  ExplicitStepperWithMassMatrix & operator=(const ExplicitStepperWithMassMatrix &) = delete;
  ExplicitStepperWithMassMatrix(ExplicitStepperWithMassMatrix &&) = default;
  ExplicitStepperWithMassMatrix & operator=(ExplicitStepperWithMassMatrix &&) = delete;
  ~ExplicitStepperWithMassMatrix() = default;

  ExplicitStepperWithMassMatrix(ode::ForwardEuler,
				SystemType && systemObj)
    : name_(StepScheme::ForwardEuler),
      order_(1),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstances_{systemObj.createRightHandSide()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrix(ode::RungeKutta4,
				SystemType && systemObj)
    : name_(StepScheme::RungeKutta4),
      order_(4),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstances_{systemObj.createRightHandSide(),
		    systemObj.createRightHandSide(),
		    systemObj.createRightHandSide(),
		    systemObj.createRightHandSide()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrix(ode::AdamsBashforth2,
				SystemType && systemObj)
    : name_(StepScheme::AdamsBashforth2),
      order_(2),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstances_{systemObj.createRightHandSide(),
                    systemObj.createRightHandSide()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrix(ode::SSPRungeKutta3,
				SystemType && systemObj)
    : name_(StepScheme::SSPRungeKutta3),
      order_(3),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstances_{systemObj.createRightHandSide()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

public:
  stepper_order_type order() const{
    return order_;
  }

  template<class LinearSolverType>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartTime,
		  ::pressio::ode::StepCount step,
		  ::pressio::ode::StepSize<independent_variable_type> dt,
		  LinearSolverType & solver)
  {

    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImpl(ode::ForwardEuler(), odeState,
		 stepStartTime.get(), dt.get(), step.get(), solver);
    }

    else if (name_ == ode::StepScheme::RungeKutta4){
      doStepImpl(ode::RungeKutta4(), odeState,
		 stepStartTime.get(), dt.get(), step.get(), solver);
    }

  }

private:
  template<class LinearSolver>
  void doStepImpl(ode::ForwardEuler,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & dt,
		  ::pressio::ode::StepCount::value_type stepNumber,
		  LinearSolver & solver)
  {

    // M(y_n, t_n, ...) (y_n+1 - y_n) = dt * rhs(t_n, ...)
    //
    // so we can do: M x = dy * rhs
    // and then y_n+1 = y_n + x

    PRESSIOLOG_DEBUG("euler forward stepper with MM: do step");

    // eval rhs and mass matrix
    auto & rhs = rhsInstances_[0];
    systemObj_.get().rightHandSide(odeState, stepStartTime, rhs);
    systemObj_.get().massMatrix(odeState, stepStartTime, massMatrix_);

    // solve for x storing into auxiliaryState_
    // CAREFUL: we are modifying the rhs in place to account for dt here
    // so pay attention afterwards if/when using the rhs
    ::pressio::ops::scale(rhs, dt);
    solver.solve(massMatrix_, auxiliaryState_, rhs);

    // need to do: y_n+1 = y_n + x
    // since odeState already contains y_n, and auxiliaryState_ = x
    // we can do: odeState += auxiliaryState_
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(odeState, one, auxiliaryState_, one);
  }

  template<class LinearSolver>
  void doStepImpl(ode::RungeKutta4,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & dt,
		  ::pressio::ode::StepCount::value_type stepNumber,
		  LinearSolver & solver)
  {
    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & rhs1 = rhsInstances_[0];
    auto & rhs2 = rhsInstances_[1];
    auto & rhs3 = rhsInstances_[2];
    auto & rhs4 = rhsInstances_[3];

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();
    constexpr auto three  = ::pressio::utils::Constants<scalar_type>::three();
    constexpr auto six  = two * three;

    const independent_variable_type dt_half = dt / two;
    const independent_variable_type t_phalf = stepStartTime + dt_half;
    const independent_variable_type dt6 = dt / six;
    const independent_variable_type dt3 = dt / three;

    auto x = ::pressio::ops::clone(odeState);
    auto auxRhs = ::pressio::ops::clone(rhs1);

    // //
    // // stage 1:
    // //
    // // rhs1 = rhs(t_n, y_n), MM = massmatrix(t_n, y_n)
    // systemObj_.get().rightHandSide(  odeState, stepStartTime, auxRhs);
    // systemObj_.get().massMatrix(odeState, stepStartTime, massMatrix_);
    // ::pressio::deep_copy(rhs1, auxRhs);
    // // solve M (ytmp - y ) = rhs1*dt_half
    // ::pressio::ops::scale(auxRhs, dt_half);
    // solver.solve(massMatrix_, x, auxRhs);
    // // ytmp = x + y
    // ::pressio::ops::update(auxiliaryState_, zero, x, one, odeState, one);

    // // stage 2:
    // // rhs2 = rhs(t_n+dt/2, ytmp), MM = massmatrix(t_n+dt/2, ytmp)
    // systemObj_.get().rightHandSide(  auxiliaryState_, t_phalf, auxRhs);
    // systemObj_.get().massMatrix(auxiliaryState_, t_phalf, massMatrix_);
    // ::pressio::deep_copy(rhs2, auxRhs);
    // ::pressio::ops::scale(auxRhs, dt_half);
    // solver.solve(massMatrix_, x, auxRhs);
    // ::pressio::ops::update(auxiliaryState_, zero, x, one, odeState, one);

    // // stage 3:
    // // rhs3 = rhs(t_n+dt/2, ytmp), MM = massmatrix(t_n+dt/2, ytmp)
    // systemObj_.get().rightHandSide(  auxiliaryState_, t_phalf, auxRhs);
    // systemObj_.get().massMatrix(auxiliaryState_, t_phalf, massMatrix_);
    // ::pressio::deep_copy(rhs3, auxRhs);
    // ::pressio::ops::scale(auxRhs, dt);
    // solver.solve(massMatrix_, x, auxRhs);
    // ::pressio::ops::update(auxiliaryState_, zero, x, one, odeState, one);

    // // stage 4:
    // systemObj_.get().rightHandSide(  auxiliaryState_, stepStartTime + dt, rhs3);
    // systemObj_.get().massMatrix(auxiliaryState_, stepStartTime + dt, massMatrix_);

    // // y_n += dt/6 * ( rhs1 + 2*rhs2 + 2*rhs3 + rhs4 )
    // ::pressio::ops::update(odeState, one,
    // 			   rhs1, dt6, rhs2, dt3, rhs3, dt3, rhs4, dt6);
  }

};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_HPP_
