/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_residual_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class GalerkinResidualType, class ProjectionPolicyType>
class ResidualPolicy : private ProjectionPolicyType
{
  mutable ::pressio::ode::step_count_type stepTracker_ = -1;

public:
  // required
  using residual_type = GalerkinResidualType;

public:
  ResidualPolicy() = delete;
  ResidualPolicy(const ResidualPolicy &) = default;
  ResidualPolicy & operator=(const ResidualPolicy &) = delete;
  ResidualPolicy(ResidualPolicy &&) = default;
  ResidualPolicy & operator=(ResidualPolicy &&) = delete;
  ~ResidualPolicy() = default;

  template<class ...Args>
  ResidualPolicy(Args && ...args)
    : ProjectionPolicyType(std::forward<Args>(args)...){}

public:
  GalerkinResidualType create() const{
    return ProjectionPolicyType::template create<GalerkinResidualType>();
  }

  template <
    class GalerkinStateType,
    class GalerkinStencilStatesContainerType,
    class GalerkinStencilVelocitiesContainerType,
    class ScalarType
    >
  void compute(::pressio::ode::SteppersE name,
	       const GalerkinStateType & galerkinState,
	       const GalerkinStencilStatesContainerType & galerkinStencilStates,
	       GalerkinStencilVelocitiesContainerType & galerkinStencilVelocities,
	       const ScalarType & t_np1,
	       const ScalarType & dt,
	       const ::pressio::ode::step_count_type & currentStepNumber,
	       GalerkinResidualType & galerkinResidual) const
  {
    if (name == ::pressio::ode::SteppersE::BDF1){
      (*this).template compute_impl_bdf<ode::BDF1>(galerkinState, galerkinStencilStates,
						   galerkinStencilVelocities, t_np1, dt,
						   currentStepNumber, galerkinResidual);
    }
    else if (name == ::pressio::ode::SteppersE::BDF2){
      (*this).template compute_impl_bdf<ode::BDF2>(galerkinState, galerkinStencilStates,
						   galerkinStencilVelocities, t_np1, dt,
						   currentStepNumber, galerkinResidual);
    }
    else if (name == ::pressio::ode::SteppersE::CrankNicolson){
      (*this).template compute_impl_cn<ode::CrankNicolson>(galerkinState, galerkinStencilStates,
							    galerkinStencilVelocities, t_np1, dt,
							    currentStepNumber, galerkinResidual);
    }
  }

private:
  template <
    class StepperTag,
    class GalerkinStateType,
    class GalerkinStencilStatesContainerType,
    class GalerkinStencilVelocitiesContainerType,
    class ScalarType
    >
  void compute_impl_bdf(const GalerkinStateType & galerkinState,
			const GalerkinStencilStatesContainerType & galerkinStencilStates,
			GalerkinStencilVelocitiesContainerType & galerkinStencilVelocities,
			const ScalarType & t_np1,
			const ScalarType & dt,
			const ::pressio::ode::step_count_type & currentStepNumber,
			GalerkinResidualType & galerkinResidual) const
  {
    ProjectionPolicyType::compute(galerkinResidual,
				  galerkinState,
				  t_np1,
				  ::pressio::ode::nPlusOne());

    ::pressio::ode::impl::discrete_time_residual(galerkinState,
						 galerkinResidual,
						 galerkinStencilStates,
						 dt, StepperTag());
  }

  template <
    class StepperTag,
    class GalerkinStateType,
    class GalerkinStencilStatesContainerType,
    class GalerkinStencilVelocitiesContainerType,
    class ScalarType
    >
  void compute_impl_cn(const GalerkinStateType & galerkinState,
		       const GalerkinStencilStatesContainerType & galerkinStencilStates,
		       GalerkinStencilVelocitiesContainerType & galerkinStencilVelocities,
		       const ScalarType & t_np1,
		       const ScalarType & dt,
		       const ::pressio::ode::step_count_type & currentStepNumber,
		       GalerkinResidualType & galerkinResidual) const
  {

    // if the step changed, I need to compute f(y_n, t_n)
    if (stepTracker_ != currentStepNumber){
      auto & f_n = galerkinStencilVelocities(::pressio::ode::n());
      auto & galState_n = galerkinStencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      ProjectionPolicyType::compute(f_n, galState_n, tn, ::pressio::ode::n());
    }

    // I always need to compute f(y_np1, t_np1)
    auto & f_np1 = galerkinStencilVelocities(::pressio::ode::nPlusOne());
    ProjectionPolicyType::compute(f_np1, galerkinState, t_np1, ::pressio::ode::nPlusOne());

    // compute discrete time residual
    ::pressio::ode::impl::discrete_time_residual
	(galerkinState, galerkinResidual, galerkinStencilStates,
	 galerkinStencilVelocities, dt, StepperTag());
  }
};

}}}}
#endif  // ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_RESIDUAL_POLICY_HPP_
