/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_jacobian_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<
  class JacobianType,
  class FomStatesManagerType,
  class DecoderType,
  class FomSystemType
  >
class UnsteadyJacobianPolicy
{

public:
  //required
  using jacobian_type = JacobianType;

public:
  UnsteadyJacobianPolicy() = delete;
  UnsteadyJacobianPolicy(const UnsteadyJacobianPolicy &) = default;
  UnsteadyJacobianPolicy & operator=(const UnsteadyJacobianPolicy &) = delete;
  UnsteadyJacobianPolicy(UnsteadyJacobianPolicy &&) = default;
  UnsteadyJacobianPolicy & operator=(UnsteadyJacobianPolicy &&) = delete;
  ~UnsteadyJacobianPolicy() = default;

  UnsteadyJacobianPolicy(const FomSystemType & fomSystem,
			 FomStatesManagerType & fomStatesMngr,
			 DecoderType & decoder)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      decoderObj_(decoder),
      decoderJacobian_(decoder.jacobianCRef())
  {}

public:
  jacobian_type create() const{
    jacobian_type J( fomSystem_.get().createApplyJacobianResult(decoderJacobian_.get()) );
    ::pressio::ops::set_zero(J);
    return J;
  }

  template <
    class LspgStateType,
    class StencilStatesContainerType,
    class ScalarType,
    class StepType
    >
  void operator()(::pressio::ode::SteppersE name,
		  const LspgStateType & lspgState,
		  const StencilStatesContainerType & stencilStates,
		  const ScalarType & time_np1,
		  const ScalarType & dt,
		  const StepType & currentStepNumber,
		  jacobian_type & lspgJacobian) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    // fomStatesMngr_.get().reconstructCurrentFomState(romState);

    // update Jacobian of decoder (if needed)
    decoderObj_.get().updateJacobian(lspgState);

    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystem_.get().applyJacobian(fomState, decoderJacobian_.get(), time_np1, lspgJacobian);

    // for lspg we have:
    //
    //		lspgJac = decoderJac + dt*coeff*J*decoderJac
    //
    // where J is the fom jacobian and coeff depends on the scheme.
    //
    // Here, lspgJacobian already contains J*decoderJac (see applyJac above),
    // so we just need to update lspgJacobian properly.

    constexpr auto one = ::pressio::utils::Constants<ScalarType>::one();
    if (name == ::pressio::ode::SteppersE::BDF1)
    {
      const auto factor = dt*::pressio::ode::constants::bdf1<ScalarType>::c_f_;
      ::pressio::ops::update(lspgJacobian, factor, decoderJacobian_.get(), one);
    }
    else if (name == ::pressio::ode::SteppersE::BDF2)
    {
      const auto factor = dt*::pressio::ode::constants::bdf2<ScalarType>::c_f_;
      ::pressio::ops::update(lspgJacobian, factor, decoderJacobian_.get(), one);
    }
    else if (name == ::pressio::ode::SteppersE::CrankNicolson)
      {
      const auto factor = dt*::pressio::ode::constants::cranknicolson<ScalarType>::c_fnp1_;
      ::pressio::ops::update(lspgJacobian, factor, decoderJacobian_.get(), one);
    }
  }

protected:
  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<DecoderType> decoderObj_ = {};
  std::reference_wrapper<const typename DecoderType::jacobian_type> decoderJacobian_ = {};
};

}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_CONTINUOUS_TIME_API_HPP_
