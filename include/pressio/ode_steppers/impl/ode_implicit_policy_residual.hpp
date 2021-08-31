/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_bdf_policy.hpp
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

#ifndef ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_
#define ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<class SystemType, class StateType, class ResidualType>
class ResidualStandardPolicyBdf
{
public:
  // required
  using residual_type = ResidualType;

public:
  ResidualStandardPolicyBdf() = delete;
  explicit ResidualStandardPolicyBdf(SystemType && systemIn)
    : systemObj_( std::forward<SystemType>(systemIn) ){}

  ResidualStandardPolicyBdf(const ResidualStandardPolicyBdf &) = default;
  ResidualStandardPolicyBdf & operator=(const ResidualStandardPolicyBdf &) = default;
  ResidualStandardPolicyBdf(ResidualStandardPolicyBdf &&) = default;
  ResidualStandardPolicyBdf & operator=(ResidualStandardPolicyBdf &&) = default;
  ~ResidualStandardPolicyBdf() = default;

public:
  ResidualType create() const{
    ResidualType R(systemObj_.get().createVelocity());
    return R;
  }

  template <
    class OdeTag,
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class ScalarType,
    class StepType
    >
  void compute(const StateType & predictedState,
	       const StencilStatesContainerType & stencilStatesManager,
	       StencilVelocitiesContainerType & stencilVelocities,
	       const ScalarType & rhsEvaluationTime,
	       const ScalarType & dt,
	       const StepType & step,
	       ResidualType & R) const
  {
    static_assert(StencilVelocitiesContainerType::size() == 0,
		  "Residual policy for BDF should have 0 velocities");

    static_assert(
		  std::is_same<OdeTag, ::pressio::ode::BDF1>::value or
		  std::is_same<OdeTag, ::pressio::ode::BDF2>::value,
		  "Invalid tag for BDF residual policy");

    try{
      systemObj_.get().velocity(predictedState, rhsEvaluationTime, R);
      ::pressio::ode::impl::discrete_time_residual(predictedState,
						   R, stencilStatesManager,
						   dt, OdeTag());
    }
    catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
};


template<class SystemType, class StateType, class ResidualType>
class ResidualStandardPolicyCrankNicolson
{
  mutable int stepTracker_ = -1;

public:
  // required
  using residual_type = ResidualType;

  ResidualStandardPolicyCrankNicolson() = delete;

  explicit ResidualStandardPolicyCrankNicolson(SystemType && systemIn)
    : systemObj_( std::forward<SystemType>(systemIn) ){}

  ResidualStandardPolicyCrankNicolson(const ResidualStandardPolicyCrankNicolson &) = default;
  ResidualStandardPolicyCrankNicolson & operator=(const ResidualStandardPolicyCrankNicolson &) = default;
  ResidualStandardPolicyCrankNicolson(ResidualStandardPolicyCrankNicolson &&) = default;
  ResidualStandardPolicyCrankNicolson & operator=(ResidualStandardPolicyCrankNicolson &&) = default;
  ~ResidualStandardPolicyCrankNicolson() = default;

public:
  ResidualType create() const{
    ResidualType R(systemObj_.get().createVelocity());
    return R;
  }

  template <
    class OdeTag,
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType,
    class ScalarType,
    class StepType
    >
  void compute(const StateType & predictedState,
    const StencilStatesContainerType & stencilStates,
    StencilVelocitiesContainerType & stencilVelocities,
    const ScalarType & t_np1,
    const ScalarType & dt,
    const StepType & step,
    ResidualType & R) const
  {
    static_assert(StencilVelocitiesContainerType::size() == 2,
      "Residual policy for CrankNicolson should have 2 velocities");

    static_assert(
      std::is_same<OdeTag, ::pressio::ode::CrankNicolson>::value,
      "Invalid tag for BDF residual policy");

    if (stepTracker_ != step){
      auto & f_n = stencilVelocities(::pressio::ode::n());
      auto & state_n = stencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      systemObj_.get().velocity(state_n, tn, f_n);
    }

    auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
    systemObj_.get().velocity(predictedState, t_np1, f_np1);
    ::pressio::ode::impl::discrete_time_residual
      (predictedState, R, stencilStates, stencilVelocities, dt, OdeTag());

    stepTracker_ = step;
  }

private:
  ::pressio::utils::InstanceOrReferenceWrapper<SystemType> systemObj_;
};

}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // ODE_IMPLICIT_IMPL_STANDARD_POLICIES_ODE_IMPLICIT_RESIDUAL_BDF_POLICY_HPP_
