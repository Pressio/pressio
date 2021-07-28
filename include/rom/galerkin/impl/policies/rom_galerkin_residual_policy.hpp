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

template <typename galerkin_residual_type, typename projection_policy_t>
class ResidualPolicy : private projection_policy_t
{
  mutable ::pressio::ode::types::step_t stepTracker_ = -1;

public:
  ResidualPolicy() = delete;
  ResidualPolicy(const ResidualPolicy &) = default;
  ResidualPolicy & operator=(const ResidualPolicy &) = delete;
  ResidualPolicy(ResidualPolicy &&) = default;
  ResidualPolicy & operator=(ResidualPolicy &&) = delete;
  ~ResidualPolicy() = default;

  template<typename ...Args>
  ResidualPolicy(std::size_t romSize, Args && ...args)
    : projection_policy_t(std::forward<Args>(args)...),
      romSize_(romSize){}

public:
  template <typename fom_system_t>
  galerkin_residual_type create(const fom_system_t & fomObj) const
  {
    galerkin_residual_type result(romSize_);
    ::pressio::ops::set_zero(result);
    return result;
  }

  template <
    typename stepper_tag,
    typename galerkin_state_t,
    typename galerkin_stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_system_t>::value
    >
  compute(const galerkin_state_t & galerkinState,
	  const galerkin_stencil_states_t & galerkinStencilStates,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  galerkin_residual_type & galerkinResidual) const
  {
    projection_policy_t::compute(galerkinResidual, galerkinState, fomSystemObj,
				 timeAtNextStep, ::pressio::ode::nPlusOne());

    ::pressio::ode::impl::discrete_time_residual
	(galerkinState, galerkinResidual,
	 galerkinStencilStates, dt, stepper_tag());
  }

  template <
    typename stepper_tag,
    typename galerkin_state_t,
    typename galerkin_stencil_states_t,
    typename stencil_velocities_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value and
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_system_t>::value
    >
  compute(const galerkin_state_t & galerkinState,
	  const galerkin_stencil_states_t & galerkinStencilStates,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & t_np1,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  stencil_velocities_t & galerkinStencilVelocities,
	  galerkin_residual_type & galerkinResidual) const
  {

    // if the step changed, I need to compute f(y_n, t_n)
    if (stepTracker_ != currentStepNumber){
      auto & f_n = galerkinStencilVelocities(::pressio::ode::n());
      auto & galState_n = galerkinStencilStates(::pressio::ode::n());
      const auto tn = t_np1-dt;
      projection_policy_t::compute(f_n, galState_n, fomSystemObj,
				   tn, ::pressio::ode::n());
    }

    // I always need to compute f(y_np1, t_np1)
    auto & f_np1 = galerkinStencilVelocities(::pressio::ode::nPlusOne());
    projection_policy_t::compute(f_np1, galerkinState, fomSystemObj,
				 t_np1, ::pressio::ode::nPlusOne());

    // compute discrete time residual
    ::pressio::ode::impl::discrete_time_residual
	(galerkinState, galerkinResidual,
	 galerkinStencilStates, galerkinStencilVelocities,
	 dt, stepper_tag());
  }

  template <
    typename stepper_tag,
    typename galerkin_state_t,
    typename galerkin_stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_discrete_time_system<fom_system_t>::value
  >
  compute(const galerkin_state_t & galerkinState,
	  const galerkin_stencil_states_t & galerkinStencilStates,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  galerkin_residual_type & galerkinResidual) const
  {
    projection_policy_t::compute(galerkinResidual, galerkinState, fomSystemObj,
				 timeAtNextStep, dt, currentStepNumber,
				 galerkinStencilStates);
  }

private:
  const std::size_t romSize_ = {};
};

}}}}
#endif  // ROM_GALERKIN_IMPL_POLICIES_ROM_GALERKIN_RESIDUAL_POLICY_HPP_
