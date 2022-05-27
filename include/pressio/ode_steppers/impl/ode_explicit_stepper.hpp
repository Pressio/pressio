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

#ifndef ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_HPP_
#define ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_HPP_

#include <array>
#include "./ode_trivial_mass_matrix.hpp"

namespace pressio{ namespace ode{ namespace impl{

// this class is not meant for direct instantiation.
// One needs to use the public create_* functions because
// templates are handled and passed properly there.
template<
  class StateType,
  class SystemType,
  class VelocityType,
  class MassMatrixType
  >
class ExplicitStepper
{

private:
  StepScheme name_;
  const stepper_order_type order_;
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  std::vector<VelocityType> velocities_;
  StateType auxiliaryState_;
  MassMatrixType massMatrix_;

public:
  ExplicitStepper() = delete;
  ExplicitStepper(const ExplicitStepper &) = default;
  ExplicitStepper & operator=(const ExplicitStepper &) = delete;
  ExplicitStepper(ExplicitStepper &&) = default;
  ExplicitStepper & operator=(ExplicitStepper &&) = delete;
  ~ExplicitStepper() = default;

  ExplicitStepper(ode::ForwardEuler,
		  SystemType && systemObj)
    : name_(StepScheme::ForwardEuler),
      order_(1),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(createMassMatrix<MassMatrixType>()(systemObj_.get()))
  {}

  ExplicitStepper(ode::RungeKutta4,
		  SystemType && systemObj)
    : name_(StepScheme::RungeKutta4),
      order_(4),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity(),
		  systemObj.createVelocity(),
		  systemObj.createVelocity(),
		  systemObj.createVelocity()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(createMassMatrix<MassMatrixType>()(systemObj_.get()))
  {}

  ExplicitStepper(ode::AdamsBashforth2,
		  SystemType && systemObj)
    : name_(StepScheme::AdamsBashforth2),
      order_(2),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity(),
                  systemObj.createVelocity()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(createMassMatrix<MassMatrixType>()(systemObj_.get()))
  {}

  ExplicitStepper(ode::SSPRungeKutta3,
		  SystemType && systemObj)
    : name_(StepScheme::SSPRungeKutta3),
      order_(3),
      systemObj_(std::forward<SystemType>(systemObj)),
      velocities_{systemObj.createVelocity()},
      auxiliaryState_{systemObj.createState()},
      massMatrix_(createMassMatrix<MassMatrixType>()(systemObj_.get()))
  {}

public:
  stepper_order_type order() const{
    return order_;
  }

  // sfinae for when mass matrix is dummy one
  template<
    class StepCountType,
    class TimeType,
    class _MassMatrixType = MassMatrixType>
  mpl::enable_if_t< is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StateType & odeState,
	     const TimeType & currentTime,
	     const TimeType & dt,
	     const StepCountType & stepNumber)
  {
    auto dummyRhsObserver = [](const StepCountType &,
			       const TimeType &,
			       const VelocityType &) { /*no op*/ };
    (*this)(odeState, currentTime, dt, stepNumber, dummyRhsObserver);
  }

  // sfinae for when mass matrix is dummy one
  template<
    class StepCountType,
    class TimeType,
    class RhsObserverType,
    class _MassMatrixType = MassMatrixType>
  mpl::enable_if_t< is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StateType & odeState,
	     const TimeType & currentTime,
	     const TimeType & dt,
	     const StepCountType & step,
	     RhsObserverType & rhsObs)
  {
    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImplNoMM(ode::ForwardEuler(), odeState,
		     currentTime, dt, step, rhsObs);
    }

    else if (name_ == ode::StepScheme::RungeKutta4){
      doStepImplNoMM(ode::RungeKutta4(), odeState,
		     currentTime, dt, step, rhsObs);
    }

    else if (name_ == ode::StepScheme::AdamsBashforth2){
      doStepImplNoMM(ode::AdamsBashforth2(), odeState,
		     currentTime, dt, step, rhsObs);
    }

    else if (name_ == ode::StepScheme::SSPRungeKutta3){
      doStepImplNoMM(ode::SSPRungeKutta3(), odeState,
		     currentTime, dt, step, rhsObs);
    }
  }

  // sfinae for when mass matrix is meaningful
  template<
    class StepCountType,
    class TimeType,
    class LinearSolverType,
    class _MassMatrixType = MassMatrixType>
  mpl::enable_if_t< !is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StateType & odeState,
	     const TimeType & currentTime,
	     const TimeType & dt,
	     const StepCountType & stepNumber,
	     LinearSolverType & solver)
  {
    auto dummyRhsObserver = [](const StepCountType &,
			       const TimeType &,
			       const VelocityType &) { /*no op*/ };
    (*this)(odeState, currentTime, dt, stepNumber,
	    solver, dummyRhsObserver);
  }

  template<
    class StepCountType,
    class TimeType,
    class LinearSolverType,
    class RhsObserverType,
    class _MassMatrixType = MassMatrixType>
  mpl::enable_if_t< !is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StateType & odeState,
	     const TimeType & currentTime,
	     const TimeType & dt,
	     const StepCountType & step,
	     LinearSolverType & solver,
	     RhsObserverType & rhsObs)
  {
    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImplWithMM(ode::ForwardEuler(), odeState,
		       currentTime, dt, step, solver, rhsObs);
    }

    // else if (name_ == ode::StepScheme::RungeKutta4){
    //   doStepImplWithMM(ode::RungeKutta4(), odeState,
    // 		       currentTime, dt, step, solver, rhsObs);
    // }

    // else if (name_ == ode::StepScheme::AdamsBashforth2){
    //   doStepImplWithMM(ode::AdamsBashforth2(), odeState,
    // 		       currentTime, dt, step, rhsObs);
    // }

    // else if (name_ == ode::StepScheme::SSPRungeKutta3){
    //   doStepImplWithMM(ode::SSPRungeKutta3(), odeState,
    // 		       currentTime, dt, step, rhsObs);
    // }
  }

private:
  template<class StepCountType, class TimeType, class RhsObserverType>
  void doStepImplNoMM(ode::ForwardEuler,
		      StateType & odeState,
		      const TimeType & currentTime,
		      const TimeType & dt,
		      const StepCountType & stepNumber,
		      RhsObserverType & rhsObs)
  {
    PRESSIOLOG_DEBUG("euler forward stepper: do step");

    auto & rhs = velocities_[0];

    //eval and observe RHS
    systemObj_.get().velocity(odeState, currentTime, rhs);
    rhsObs(stepNumber, currentTime, rhs);

    // y = y + dt * rhs
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(odeState, one, rhs, dt);
  }

  template<
    class StepCountType, class TimeType,
    class LinearSolver, class RhsObserverType
    >
  void doStepImplWithMM(ode::ForwardEuler,
			StateType & odeState,
			const TimeType & currentTime,
			const TimeType & dt,
			const StepCountType & stepNumber,
			LinearSolver & solver,
			RhsObserverType & rhsObs)
  {
    //
    // M(t_n, ...) (y_n+1 - y_n) = dt * rhs(t_n, ...)
    //
    // so we can do: M x = dy * rhs
    // and then y_n+1 = y_n + x

    PRESSIOLOG_DEBUG("euler forward stepper: do step");

    // eval rhs and mass matrix
    auto & rhs = velocities_[0];
    systemObj_.get().velocity(odeState, currentTime, rhs);
    systemObj_.get().massMatrix(odeState, currentTime, massMatrix_);

    // observe RHS
    rhsObs(stepNumber, currentTime, rhs);

    // solve for x storing into auxiliaryState_
    // CAREFUL: we are changing rhs here so pay attention afterwards
    ::pressio::ops::scale(rhs, dt);
    solver.solve(massMatrix_, auxiliaryState_, rhs);

    // need to do: y_n+1 = y_n + x
    // since y_n+1 already contains y_n, and auxiliaryState_ = x
    // we can do: y_n+1 += auxiliaryState_
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(odeState, one, auxiliaryState_, one);
  }

  template<class StepCountType, class TimeType, class RhsObserverType>
  void doStepImplNoMM(ode::AdamsBashforth2,
		  StateType & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const StepCountType & stepNumber,
		  RhsObserverType & rhsObs)
  {
    PRESSIOLOG_DEBUG("adams-bashforth2 stepper: do step");

    // y_n+1 = y_n + dt*[ (3/2)*f(y_n, t_n) - (1/2)*f(y_n-1, t_n-1) ]

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    const auto cfn   = ::pressio::utils::Constants<scalar_type>::threeOvTwo()*dt;
    const auto cfnm1 = ::pressio::utils::Constants<scalar_type>::negOneHalf()*dt;

    if (stepNumber==1){
      // use Euler forward or we could use something else here maybe RK4
      auto & rhs = velocities_[0];
      systemObj_.get().velocity(odeState, currentTime, rhs);
      rhsObs(stepNumber, currentTime, rhs);

      using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
      ::pressio::ops::update(odeState, one, rhs, dt);
    }
    else{
      auto & fn   = velocities_[0];
      auto & fnm1 = velocities_[1];
      // fn -> fnm1
      ::pressio::ops::deep_copy(fnm1, fn);

      systemObj_.get().velocity(odeState, currentTime, fn);
      rhsObs(stepNumber, currentTime, fn);

      using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
      constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
      ::pressio::ops::update(odeState, one, fn, cfn, fnm1, cfnm1);
    }
  }

  template<class StepCountType, class TimeType, class RhsObserverType>
  void doStepImplNoMM(ode::SSPRungeKutta3,
		  StateType & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const StepCountType & stepNumber,
		  RhsObserverType & rhsObs)
  {
    PRESSIOLOG_DEBUG("ssprk3 stepper: do step");

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one   = ::pressio::utils::Constants<scalar_type>::one();
    constexpr auto two   = ::pressio::utils::Constants<scalar_type>::two();
    constexpr auto three = ::pressio::utils::Constants<scalar_type>::three();
    constexpr auto four  = ::pressio::utils::Constants<scalar_type>::four();
    constexpr auto oneOvTwo = one/two;
    constexpr auto oneOvThree = one/three;
    constexpr auto twoOvThree = two/three;
    constexpr auto threeOvFour = three/four;
    constexpr auto fourInv = one/four;

    auto & rhs0 = velocities_[0];

    // see e.g. https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html

    // rhs(u_n, t_n)
    systemObj_.get().velocity(odeState, currentTime, rhs0);
    // u_1 = u_n + dt * rhs(u_n, t_n)
    ::pressio::ops::update(auxiliaryState_, zero,
                           odeState,        one,
                           rhs0,            dt);
    rhsObs(stepNumber, currentTime, rhs0);

    // rhs(u_1, t_n+dt)
    systemObj_.get().velocity(auxiliaryState_, currentTime+dt, rhs0);
    // u_2 = 3/4*u_n + 1/4*u_1 + 1/4*dt*rhs(u_1, t_n+dt)
    ::pressio::ops::update(auxiliaryState_, fourInv,
			   odeState,        threeOvFour,
                           rhs0,            fourInv*dt);

    // rhs(u_2, t_n + 0.5*dt)
    systemObj_.get().velocity(auxiliaryState_, currentTime + oneOvTwo*dt, rhs0);
    // u_n+1 = 1/3*u_n + 2/3*u_2 + 2/3*dt*rhs(u_2, t_n+0.5*dt)
    ::pressio::ops::update(odeState,        oneOvThree,
			   auxiliaryState_, twoOvThree,
                           rhs0,            twoOvThree*dt);
  }

  template<class StepCountType, class TimeType, class RhsObserverType>
  void doStepImplNoMM(ode::RungeKutta4,
		  StateType & odeState,
		  const TimeType & currentTime,
		  const TimeType & dt,
		  const StepCountType & stepNumber,
		  RhsObserverType & rhsObs)
  {
    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & rhs0 = velocities_[0];
    auto & rhs1 = velocities_[1];
    auto & rhs2 = velocities_[2];
    auto & rhs3 = velocities_[3];

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto two  = ::pressio::utils::Constants<scalar_type>::two();
    constexpr auto three  = ::pressio::utils::Constants<scalar_type>::three();
    constexpr auto six  = two * three;

    const TimeType dt_half = dt / two;
    const TimeType t_phalf = currentTime + dt_half;
    const TimeType dt6 = dt / six;
    const TimeType dt3 = dt / three;

    // stage 1: ytmp = y + rhs0*dt_half;
    systemObj_.get().velocity(odeState, currentTime, rhs0);
    rhsObs(stepNumber, currentTime, rhs0);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs0, dt_half);

    // stage 2: ytmp = y + rhs1*dt_half;
    systemObj_.get().velocity(auxiliaryState_, t_phalf, rhs1);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs1, dt_half);

    // stage 3: ytmp = y + rhs2*dt;
    systemObj_.get().velocity(auxiliaryState_, t_phalf, rhs2);
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    systemObj_.get().velocity(auxiliaryState_, currentTime + dt, rhs3);
    this->rk4_stage_update_impl(odeState, rhs0, rhs1, rhs2, rhs3, dt6, dt3);
  }

  template<class rhs_t, class TimeType>
  void rk4_stage_update_impl(StateType & yIn,
			     const StateType & stateIn,
			     const rhs_t & rhsIn,
			     TimeType dtValue)
  {
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto zero  = ::pressio::utils::Constants<scalar_type>::zero();
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(yIn, zero, stateIn, one, rhsIn, dtValue);
  }

  template<class rhs_t, class TimeType>
  void rk4_stage_update_impl(StateType & stateIn,
			     const rhs_t & rhsIn0,
			     const rhs_t & rhsIn1,
			     const rhs_t & rhsIn2,
			     const rhs_t & rhsIn3,
			     TimeType dt6,
			     TimeType dt3)
  {
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one  = ::pressio::utils::Constants<scalar_type>::one();
    ::pressio::ops::update(stateIn, one,
			   rhsIn0, dt6,
			   rhsIn1, dt3,
			   rhsIn2, dt3,
			   rhsIn3, dt6);
  }

};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_STEPPER_HPP_
