/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_velocity_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  class ScalarType,
  class GalerkinStateType,
  class GalerkinVelocityType,
  class ProjectionPolicyType
  >
class VelocityOnlySystemUsing
{
public:
  // required
  using scalar_type   = ScalarType;
  using state_type    = GalerkinStateType;
  using velocity_type = GalerkinVelocityType;

  VelocityOnlySystemUsing() = delete;
  VelocityOnlySystemUsing(const VelocityOnlySystemUsing &) = default;
  VelocityOnlySystemUsing & operator=(const VelocityOnlySystemUsing &) = delete;
  VelocityOnlySystemUsing(VelocityOnlySystemUsing &&) = default;
  VelocityOnlySystemUsing & operator=(VelocityOnlySystemUsing &&) = delete;
  ~VelocityOnlySystemUsing() = default;

  template<class ...Args>
  VelocityOnlySystemUsing(const GalerkinStateType & romState,
			  Args && ...args)
    : pp_(::pressio::ops::extent(romState,0), std::forward<Args>(args)...)
  {}

public:
  GalerkinVelocityType createVelocity() const
  {
    return pp_.template create<GalerkinVelocityType>();
  }

  void velocity(const GalerkinStateType & galerkinState,
		const ScalarType & timeForRhsEvaluation,
		GalerkinVelocityType & galerkinRhs) const
  {
    pp_.compute(galerkinRhs, galerkinState,
		timeForRhsEvaluation, ::pressio::ode::n());
  }
private:
    ProjectionPolicyType pp_;
};

template <
  class ScalarType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class ProjectionPolicyType1,
  class ProjectionPolicyType2
  >
struct DiscreteTimeReducedSystem
{
  using scalar_type = ScalarType;
  using state_type  = GalerkinStateType;
  using discrete_time_residual_type = GalerkinResidualType;
  using discrete_time_jacobian_type = GalerkinJacobianType;

  DiscreteTimeReducedSystem() = delete;
  DiscreteTimeReducedSystem(const DiscreteTimeReducedSystem &) = default;
  DiscreteTimeReducedSystem & operator=(const DiscreteTimeReducedSystem &) = delete;
  DiscreteTimeReducedSystem(DiscreteTimeReducedSystem &&) = default;
  DiscreteTimeReducedSystem & operator=(DiscreteTimeReducedSystem &&) = delete;
  ~DiscreteTimeReducedSystem() = default;

  template<class T1, class T2, class T3, class T4, class T5>
  DiscreteTimeReducedSystem(const T1 & romStateIn,
			    const T2 & projector,
			    const T3 & fomObj,
			    const T4 & decoder,
			    T5 & fomStatesMngr)
    : p1(::pressio::ops::extent(romStateIn,0),
	 projector, fomObj, fomStatesMngr),
      // this is on purpose because for Galerkin we have square jacobian
      p2(::pressio::ops::extent(romStateIn,0),
	 ::pressio::ops::extent(romStateIn,0),
	 projector, fomObj, fomStatesMngr, decoder)
  {}

  template<class T1, class T2, class T3, class T4, class T5, class T6>
  DiscreteTimeReducedSystem(const T1 & romStateIn,
			    const T2 & projector,
			    const T3 & masker,
			    const T4 & fomObj,
			    const T5 & decoder,
			    T6 & fomStatesMngr)
    : p1(::pressio::ops::extent(romStateIn,0),
	 projector, masker, fomObj, fomStatesMngr),
      // this is on purpose because for Galerkin we have square jacobian
      p2(::pressio::ops::extent(romStateIn,0),
	 ::pressio::ops::extent(romStateIn,0),
	 projector, masker, fomObj, fomStatesMngr, decoder)
  {}

  discrete_time_residual_type createDiscreteTimeResidual() const{
    return p1.template create<discrete_time_residual_type>();
  }

  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
    return p2.template create<discrete_time_jacobian_type>();
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType step_count,
                            scalar_type time_np1,
                            scalar_type dt,
                            discrete_time_residual_type & R,
                            const state_type & gal_state_np1) const
  {
    p1.compute(R, gal_state_np1, time_np1, dt, step_count);
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType step_count,
                            scalar_type time_np1,
                            scalar_type dt,
                            discrete_time_residual_type & R,
                            const state_type & gal_state_np1,
			    const state_type & gal_state_n) const
  {
    p1.compute(R, gal_state_np1, time_np1, dt,
    	       step_count, gal_state_n);
  }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType step_count,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_residual_type & R,
			    const state_type & gal_state_np1,
			    const state_type & gal_state_n,
			    const state_type & gal_state_nm1) const
  {
    p1.compute(R, gal_state_np1, time_np1, dt, step_count,
     	       gal_state_n, gal_state_nm1);
  }

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType step_count,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    const state_type & gal_state_np1) const
  {
    p2.compute(J, gal_state_np1, time_np1, dt, step_count);
  }

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType step_count,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    const state_type & gal_state_np1,
			    const state_type & gal_state_n) const
  {
    p2.compute(J, gal_state_np1, time_np1, dt, step_count,gal_state_n);
  }

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType step_count,
			    scalar_type time_np1,
			    scalar_type dt,
			    discrete_time_jacobian_type & J,
			    const state_type & gal_state_np1,
			    const state_type & gal_state_n,
			    const state_type & gal_state_nm1) const
  {
    p2.compute(J, gal_state_np1, time_np1, dt,
	       step_count, gal_state_n, gal_state_nm1);
  }

private:
  ProjectionPolicyType1 p1;
  ProjectionPolicyType2 p2;

};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_VELOCITY_POLICY_HPP_
