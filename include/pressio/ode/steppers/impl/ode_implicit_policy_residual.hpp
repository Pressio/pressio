/*
//@HEADER
// ************************************************************************
//
// ode_implicit_policy_residual.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_HPP_
#define ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_HPP_

#include "./ode_trivial_mass_matrix.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<
  class SystemType,
  class TimeType,
  class StateType,
  class ResidualType,
  class MassMatrixType = NoOpMassMatrix
  >
class ResidualStandardPolicy
{
public:
  // required
  using time_type     = TimeType;
  using state_type    = StateType;
  using residual_type = ResidualType;

public:
  ResidualStandardPolicy() = delete;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  explicit ResidualStandardPolicy(SystemType systemIn)
    : systemObj_(systemIn),
      massMatrix_(createMassMatrix<MassMatrixType>()(systemIn)){}
#else
  explicit ResidualStandardPolicy(SystemType && systemIn)
    : systemObj_( std::forward<SystemType>(systemIn) ),
      massMatrix_(createMassMatrix<MassMatrixType>()(systemIn)){}
#endif

  ResidualStandardPolicy(const ResidualStandardPolicy &) = default;
  ResidualStandardPolicy & operator=(const ResidualStandardPolicy &) = default;
  ResidualStandardPolicy(ResidualStandardPolicy &&) = default;
  ResidualStandardPolicy & operator=(ResidualStandardPolicy &&) = default;
  ~ResidualStandardPolicy() = default;

public:
  StateType createState() const{
    StateType result(systemObj_.get().createState());
    return result;
  }

  ResidualType create() const{
    ResidualType R(systemObj_.get().createVelocity());
    return R;
  }

  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class StepType,
    class _MassMatrixType = MassMatrixType>
  mpl::enable_if_t< is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StepScheme name,
	     const StateType & predictedState,
	     const StencilStatesContainerType & stencilStatesManager,
	     StencilVelocitiesContainerType & stencilVelocities,
	     const TimeType & rhsEvaluationTime,
	     const TimeType & dt,
	     const StepType & step,
	     ResidualType & R) const
  {

    if (name == StepScheme::BDF1){
      (*this).template compute_impl_bdf_no_mm<ode::BDF1>(predictedState, stencilStatesManager,
							 stencilVelocities, rhsEvaluationTime,
							 dt, step, R);
    }
    else if (name == StepScheme::BDF2){
      (*this).template compute_impl_bdf_no_mm<ode::BDF2>(predictedState, stencilStatesManager,
							 stencilVelocities, rhsEvaluationTime,
							 dt, step, R);
    }
    else if (name == StepScheme::CrankNicolson){
      this->compute_impl_cn_no_mm(predictedState, stencilStatesManager,
				  stencilVelocities, rhsEvaluationTime, dt, step, R);
    }
  }

  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class StepType,
    class _MassMatrixType = MassMatrixType
    >
  mpl::enable_if_t< !is_trivial_mass_matrix<_MassMatrixType>::value >
  operator()(StepScheme name,
	     const StateType & predictedState,
	     const StencilStatesContainerType & stencilStatesManager,
	     StencilVelocitiesContainerType & stencilVelocities,
	     const TimeType & rhsEvaluationTime,
	     const TimeType & dt,
	     const StepType & step,
	     ResidualType & R) const
  {

    if (name == StepScheme::BDF1){
      (*this).template compute_impl_bdf_with_mm<ode::BDF1>(predictedState, stencilStatesManager,
							 stencilVelocities, rhsEvaluationTime,
							 dt, step, R);
    }
    // else if (name == StepScheme::BDF2){
    //   (*this).template compute_impl_bdf_with_mm<ode::BDF2>(predictedState, stencilStatesManager,
    // 							 stencilVelocities, rhsEvaluationTime,
    // 							 dt, step, R);
    // }
    // else if (name == StepScheme::CrankNicolson){
    //   this->compute_impl_cn_with_mm(predictedState, stencilStatesManager,
    // 				  stencilVelocities, rhsEvaluationTime, dt, step, R);
    // }
  }

private:
  //
  // BDF NO mass matrix
  //
  template <
  class OdeTag,
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepType
  >
  void compute_impl_bdf_no_mm(const StateType & predictedState,
			      const StencilStatesContainerType & stencilStatesManager,
			      StencilVelocitiesContainerType & stencilVelocities,
			      const TimeType & rhsEvaluationTime,
			      const TimeType & dt,
			      const StepType & step,
			      ResidualType & R) const
  {

    try{
      systemObj_.get().velocity(predictedState, rhsEvaluationTime, R);
      ::pressio::ode::impl::discrete_time_residual(predictedState,
						   R, stencilStatesManager,
						   dt, OdeTag());
      stepTracker_ = step;
    }
    catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  //
  // BDF WITH mass matrix
  //
  template <
  class OdeTag,
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepType
  >
  void compute_impl_bdf_with_mm(const StateType & predictedState,
				const StencilStatesContainerType & stencilStatesManager,
				StencilVelocitiesContainerType & stencilVelocities,
				const TimeType & rhsEvaluationTime,
				const TimeType & dt,
				const StepType & step,
				ResidualType & R) const
  {

    // try{
    // systemObj_.get().velocity(predictedState, rhsEvaluationTime, R);
    // systemObj_.get().massMatrix(predictedState, rhsEvaluationTime, massMatrix_);
    // ::pressio::ode::impl::discrete_time_residual(predictedState,
    // 						   R, stencilStatesManager,
    // 						   dt, OdeTag());
    //   stepTracker_ = step;
    // }
    // catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
    //   throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    // }
  }

  //
  // CN
  //
  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class StepType
    >
  void compute_impl_cn_no_mm(const StateType & predictedState,
			     const StencilStatesContainerType & stencilStates,
			     StencilVelocitiesContainerType & stencilVelocities,
			     const TimeType & t_np1,
			     const TimeType & dt,
			     const StepType & step,
			     ResidualType & R) const
  {

    if (stepTracker_ != step){
      auto & f_n = stencilVelocities(::pressio::ode::n());
      auto & state_n = stencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      systemObj_.get().velocity(state_n, tn, f_n);
    }

    auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
    systemObj_.get().velocity(predictedState, t_np1, f_np1);
    ::pressio::ode::impl::discrete_time_residual
      (predictedState, R, stencilStates, stencilVelocities, dt,
       ode::CrankNicolson());

    stepTracker_ = step;
  }

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  mutable int32_t stepTracker_ = -1;
  MassMatrixType massMatrix_;
};

}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_HPP_
