

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

#ifndef ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_HPP_
#define ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class SystemType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class JacobianType
  >
class ResidualJacobianStandardPolicy
{

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = StateType;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

public:
  ResidualJacobianStandardPolicy() = delete;

  explicit ResidualJacobianStandardPolicy(SystemType && systemIn)
    : systemObj_( std::forward<SystemType>(systemIn) ){}

  ResidualJacobianStandardPolicy(const ResidualJacobianStandardPolicy &) = default;
  ResidualJacobianStandardPolicy & operator=(const ResidualJacobianStandardPolicy &) = default;
  // ResidualJacobianStandardPolicy(ResidualJacobianStandardPolicy &&) = default;
  // ResidualJacobianStandardPolicy & operator=(ResidualJacobianStandardPolicy &&) = default;
  ~ResidualJacobianStandardPolicy() = default;

public:
  StateType createState() const{
    StateType result(systemObj_.get().createState());
    return result;
  }

  ResidualType createResidual() const{
    ResidualType R(systemObj_.get().createRightHandSide());
    return R;
  }

  JacobianType createJacobian() const{
    JacobianType JJ(systemObj_.get().createJacobian());
    return JJ;
  }

  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType>
  void operator()(StepScheme name,
		  const StateType & predictedState,
		  const StencilStatesContainerType & stencilStatesManager,
		  StencilVelocitiesContainerType & stencilVelocities,
		  const ::pressio::ode::StepEndAt<IndVarType> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  ResidualType & R,
		  JacobianType & J,
		  bool computeJacobian) const
  {

    if (name == StepScheme::BDF1){
      (*this).template compute_impl_bdf
	<ode::BDF1>(predictedState, stencilStatesManager,
		    stencilVelocities, rhsEvaluationTime.get(),
		    dt.get(), step.get(), R, J, computeJacobian);
    }

    else if (name == StepScheme::BDF2){
      (*this).template compute_impl_bdf
	<ode::BDF2>(predictedState, stencilStatesManager,
		    stencilVelocities, rhsEvaluationTime.get(),
		    dt.get(), step.get(), R, J, computeJacobian);
    }

    else if (name == StepScheme::CrankNicolson){
      this->compute_impl_cn(predictedState, stencilStatesManager,
			    stencilVelocities, rhsEvaluationTime.get(),
			    dt.get(), step.get(), R, J, computeJacobian);
    }
  }

private:
  //
  // BDF
  //
  template <
  class OdeTag,
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepType
  >
  void compute_impl_bdf(const StateType & predictedState,
			const StencilStatesContainerType & stencilStatesManager,
			StencilVelocitiesContainerType & /*unused*/,
			const IndVarType & rhsEvaluationTime,
			const IndVarType & dt,
			const StepType & step,
			ResidualType & R,
			JacobianType & J,
			bool computeJacobian) const
  {

    try{
      systemObj_.get()(predictedState, rhsEvaluationTime, R, J, computeJacobian);
      ::pressio::ode::impl::discrete_residual(OdeTag(), predictedState,
					      R, stencilStatesManager, dt);

      if (computeJacobian){
	::pressio::ode::impl::discrete_jacobian(OdeTag(), J, dt);
      }

      stepTracker_ = step;
    }
    catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  //
  // CN
  //
  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class StepType
    >
  void compute_impl_cn(const StateType & predictedState,
		       const StencilStatesContainerType & stencilStates,
		       StencilVelocitiesContainerType & stencilVelocities,
		       const IndVarType & t_np1,
		       const IndVarType & dt,
		       const StepType & step,
		       ResidualType & R,
		       JacobianType & J,
		       bool computeJacobian) const
  {

    if (stepTracker_ != step){
      auto & f_n = stencilVelocities(::pressio::ode::n());
      auto & state_n = stencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      systemObj_.get()(state_n, tn, f_n, J, false);
    }

    auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
    systemObj_.get()(predictedState, t_np1, f_np1, J, computeJacobian);
    ::pressio::ode::impl::discrete_residual
	(ode::CrankNicolson(), predictedState, R, stencilStates,
	 stencilVelocities, dt);

    if (computeJacobian){
      //systemObj_.get().jacobian(predictedState, t_np1, J);
      ::pressio::ode::impl::discrete_jacobian(ode::CrankNicolson(), J, dt);
    }

    stepTracker_ = step;
  }

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
  mutable int32_t stepTracker_ = -1;
};

}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_HPP_
