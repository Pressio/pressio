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

#include <vector>

namespace pressio{ namespace ode{ namespace impl{

// this class is NOT meant for direct instantiation.
// One needs to use the public create_* functions because
// templates are handled and passed properly there.

template<
  class StateType,
  class IndVarType,
  class SystemType,
  class RightHandSideType
  >
class ExplicitStepperWithMassMatrixImpl
{

public:
  using independent_variable_type  = IndVarType;
  using state_type  = StateType;

private:
  StepScheme name_;
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  RightHandSideType rhsInstance_;

  // xInstances is a container of instances of states
  // that are used in the solve M x = b
  std::vector<StateType> xInstances_;

  typename mpl::remove_cvref_t<SystemType>::mass_matrix_type massMatrix_;

public:
  ExplicitStepperWithMassMatrixImpl() = delete;
  ExplicitStepperWithMassMatrixImpl(const ExplicitStepperWithMassMatrixImpl &) = default;
  ExplicitStepperWithMassMatrixImpl & operator=(const ExplicitStepperWithMassMatrixImpl &) = delete;
  ~ExplicitStepperWithMassMatrixImpl() = default;

  ExplicitStepperWithMassMatrixImpl(ode::ForwardEuler /*tag*/,
				    SystemType && systemObj)
    : name_(StepScheme::ForwardEuler),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstance_{systemObj.createRightHandSide()},
      xInstances_{systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrixImpl(ode::RungeKutta4  /*tag*/,
				    SystemType && systemObj)
    : name_(StepScheme::RungeKutta4),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstance_{systemObj.createRightHandSide()},
      xInstances_{systemObj.createState(),
		  systemObj.createState(),
		  systemObj.createState(),
		  systemObj.createState(),
		  systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrixImpl(ode::AdamsBashforth2 /*tag*/,
				    SystemType && systemObj)
    : name_(StepScheme::AdamsBashforth2),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstance_{systemObj.createRightHandSide()},
      xInstances_{systemObj.createState(),
                  systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

  ExplicitStepperWithMassMatrixImpl(ode::SSPRungeKutta3 /*tag*/,
				    SystemType && systemObj)
    : name_(StepScheme::SSPRungeKutta3),
      systemObj_(std::forward<SystemType>(systemObj)),
      rhsInstance_{systemObj.createRightHandSide()},
      xInstances_{systemObj.createState(),
                  systemObj.createState()},
      massMatrix_(systemObj.createMassMatrix())
  {}

public:
  template<class LinearSolverType>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount step,
		  ::pressio::ode::StepSize<independent_variable_type> stepSize,
		  LinearSolverType & solver)
  {
    auto dummyRhsObserver = [](::pressio::ode::StepCount /*unused*/,
			       ::pressio::ode::IntermediateStepCount /*unused*/,
			       const independent_variable_type & /*unused*/,
			       const RightHandSideType & /*unused*/)
    {
      /*no op*/
    };

    (*this)(odeState, stepStartVal, step, stepSize, solver, dummyRhsObserver);
  }

  template<class LinearSolverType, class RhsObserverType>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount step,
		  ::pressio::ode::StepSize<independent_variable_type> stepSize,
		  LinearSolverType & solver,
		  RhsObserverType & rhsObserver)
  {

    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImpl(ode::ForwardEuler(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, solver, rhsObserver);
    }

    else if (name_ == ode::StepScheme::AdamsBashforth2){
      doStepImpl(ode::AdamsBashforth2(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, solver, rhsObserver);
    }

    else if (name_ == ode::StepScheme::RungeKutta4){
      doStepImpl(ode::RungeKutta4(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, solver, rhsObserver);
    }

    else if (name_ == ode::StepScheme::SSPRungeKutta3){
      doStepImpl(ode::SSPRungeKutta3(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, solver, rhsObserver);
    }
  }

private:

  template<class LinearSolver, class RhsObserverType>
  void doStepImpl(ode::ForwardEuler,
		  StateType & odeState,
		  const independent_variable_type & stepStartVal,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  LinearSolver & solver,
		  RhsObserverType & rhsObserver)
  {

    // M(y_n, t_n, ...) (y_n+1 - y_n) = stepSize * rhs(t_n, ...)
    //
    // so we do: M x = f, and then y_n+1 = y_n + x*stepSize

    PRESSIOLOG_DEBUG("euler forward stepper with MM: do step");

    auto & fn = rhsInstance_;
    auto & x  = xInstances_[0];

    systemObj_.get()(odeState, stepStartVal, fn, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartVal, fn);
    solver.solve(massMatrix_, x, fn);

    // need to do: y_n+1 = y_n + stepSize*x
    // odeState already contains y_n
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(odeState, one, x, stepSize);
  }

  template<class LinearSolver, class RhsObserverType>
  void doStepImpl(ode::AdamsBashforth2,
		  StateType & odeState,
		  const independent_variable_type & stepStartVal,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  LinearSolver & solver,
		  RhsObserverType & rhsObserver)
  {

    /*
      y_n+1 = y_n + (3/2)*h*M_n^-1*f_n - (1/2)*h*M_n-1^-1*f_n-1
     */

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();

    auto & fn = rhsInstance_;

    if (stepNumber.get()==1){
      // start up with Euler forward

      auto & x  = xInstances_[0];
      systemObj_.get()(odeState, stepStartVal, fn, massMatrix_);
      rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartVal, fn);
      solver.solve(massMatrix_, x, fn);

      // now compute new state y_n+1 = y_n + dt * x
      ::pressio::ops::update(odeState, one, x, stepSize);
    }
    else{
      auto & xn   = xInstances_[0];
      auto & xnm1 = xInstances_[1];
      ::pressio::ops::deep_copy(xnm1, xn);

      systemObj_.get()(odeState, stepStartVal, fn, massMatrix_);
      rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartVal, fn);
      solver.solve(massMatrix_, xn, fn);

      const auto cfn   = ::pressio::utils::Constants<scalar_type>::threeOvTwo()*stepSize;
      const auto cfnm1 = ::pressio::utils::Constants<scalar_type>::negOneHalf()*stepSize;
      ::pressio::ops::update(odeState, one, xn, cfn, xnm1, cfnm1);
    }

  }

  template<class LinearSolver, class RhsObserverType>
  void doStepImpl(ode::SSPRungeKutta3,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  LinearSolver & solver,
		  RhsObserverType & rhsObserver)
  {
    PRESSIOLOG_DEBUG("ssprk3 stepper: do step");

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one   = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two   = ::pressio::utils::Constants<scalar_type>::two();
    constexpr auto three = ::pressio::utils::Constants<scalar_type>::three();
    constexpr auto four  = ::pressio::utils::Constants<scalar_type>::four();
    constexpr auto oneOvThree = one/three;
    constexpr auto twoOvThree = two/three;
    constexpr auto threeOvFour = three/four;
    constexpr auto fourInv = one/four;

    auto & rhs = rhsInstance_;
    auto & x = xInstances_[0];
    auto & auxState = xInstances_[1];

    // see e.g. https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html

    const scalar_type stepSize_half{stepSize/two};
    const independent_variable_type t_phalf{stepStartTime + stepSize_half};
    const independent_variable_type t_next{stepStartTime + stepSize};

    // rhs(u_n, t_n)
    systemObj_.get()(odeState, stepStartTime, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs);
    solver.solve(massMatrix_, x, rhs);
    // u_1 = u_n + stepSize * x
    ::pressio::ops::update(auxState, zero, odeState, one, x, stepSize);

    // rhs(u_1, t_n+stepSize)
    systemObj_.get()(auxState, t_next, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(1), t_next, rhs);
    solver.solve(massMatrix_, x, rhs);
    // u_2 = 3/4*u_n + 1/4*u_1 + 1/4*stepSize*x
    ::pressio::ops::update(auxState, fourInv, odeState, threeOvFour, x, fourInv*stepSize);

    // rhs(u_2, t_n + 0.5*stepSize)
    systemObj_.get()(auxState, t_phalf, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(2), t_phalf, rhs);
    solver.solve(massMatrix_, x, rhs);
    // u_n+1 = 1/3*u_n + 2/3*u_2 + 2/3*stepSize*rhs(u_2, t_n+0.5*stepSize)
    ::pressio::ops::update(odeState, oneOvThree, auxState, twoOvThree, x, twoOvThree*stepSize);
  }

  template<class LinearSolver, class RhsObserverType>
  void doStepImpl(ode::RungeKutta4,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  LinearSolver & solver,
		  RhsObserverType & rhsObserver)
  {

    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & rhs = rhsInstance_;
    auto & x1 = xInstances_[0];
    auto & x2 = xInstances_[1];
    auto & x3 = xInstances_[2];
    auto & x4 = xInstances_[3];
    auto & auxState = xInstances_[4];

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();
    constexpr auto three  = ::pressio::utils::Constants<scalar_type>::three();
    constexpr auto six  = two * three;

    const scalar_type stepSize_half{stepSize / two};
    const independent_variable_type t_phalf{stepStartTime + stepSize_half};
    const independent_variable_type t_next{stepStartTime + stepSize};
    const scalar_type stepSize6{stepSize / six};
    const scalar_type stepSize3{stepSize / three};

    // stage 1:
    // rhs1 = rhs(y_n, t_n)
    systemObj_.get()(odeState, stepStartTime, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs);
    solver.solve(massMatrix_, x1, rhs);

    // stage 2:
    // ytmp = y + rhs1*stepSize_half;
    this->rk4_stage_update_impl(auxState, odeState, x1, stepSize_half);
    // rhs2 = rhs(y_tmp, t_n+stepSize/2)
    systemObj_.get()(auxState, t_phalf, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(1), t_phalf, rhs);
    solver.solve(massMatrix_, x2, rhs);

    // stage 3:
    // ytmp = y + rhs2*stepSize_half;
    this->rk4_stage_update_impl(auxState, odeState, x2, stepSize_half);
    // rhs3 = rhs(y_tmp)
    systemObj_.get()(auxState, t_phalf, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(2), t_phalf, rhs);
    solver.solve(massMatrix_, x3, rhs);

    // stage 4:
    // ytmp = y + rhs3*stepSize;
    this->rk4_stage_update_impl(auxState, odeState, x3, stepSize);
    // rhs3 = rhs(y_tmp)
    systemObj_.get()(auxState, t_next, rhs, massMatrix_);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(3), t_next, rhs);
    solver.solve(massMatrix_, x4, rhs);

    ::pressio::ops::update(odeState, one,
			   x1, stepSize6, x2, stepSize3,
			   x3, stepSize3, x4, stepSize6);
  }

  template<class Xt, class FactorType>
  void rk4_stage_update_impl(StateType & yIn,
			     const StateType & stateIn,
			     const Xt & xIn,
			     const FactorType & rhsFactor)
  {
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(yIn, zero, stateIn, one, xIn, rhsFactor);
  }

};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_WITH_MASS_MATRIX_HPP_
