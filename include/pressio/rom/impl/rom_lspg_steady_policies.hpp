/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_policies.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_STEADY_POLICIES_HPP_
#define ROM_IMPL_ROM_LSPG_STEADY_POLICIES_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class FomStatesManagerType, class ResidualType, class FomSystemType>
class SteadyResidualPolicy
{
public:
  using residual_type = ResidualType;

public:
  SteadyResidualPolicy() = delete;
  SteadyResidualPolicy(const SteadyResidualPolicy &) = default;
  SteadyResidualPolicy & operator=(const SteadyResidualPolicy &) = delete;
  SteadyResidualPolicy(SteadyResidualPolicy &&) = default;
  SteadyResidualPolicy & operator=(SteadyResidualPolicy &&) = delete;
  ~SteadyResidualPolicy() = default;

  SteadyResidualPolicy(const FomSystemType & fomSystem,
		       FomStatesManagerType & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr), fomSystem_(fomSystem){}

public:
  residual_type create() const
  {
    return residual_type(fomSystem_.get().createResidual());
  }

  template <class LspgStateType>
  void compute(const LspgStateType & lspgState,
	       residual_type & romResidual) const
  {
    fomStatesMngr_.get().reconstructCurrentFomState(lspgState);
    const auto & currFomState = fomStatesMngr_.get().currentFomState();
    fomSystem_.get().residual(currFomState, romResidual);
  }

protected:
  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
};

template<
  class FomStatesManagerType,
  class JacobianType,
  class DecoderType,
  class FomSystemType
  >
class SteadyJacobianPolicy
{
public:
  using jacobian_type = JacobianType;

public:
  SteadyJacobianPolicy() = delete;
  SteadyJacobianPolicy(const SteadyJacobianPolicy &) = default;
  SteadyJacobianPolicy & operator=(const SteadyJacobianPolicy &) = delete;
  SteadyJacobianPolicy(SteadyJacobianPolicy &&) = default;
  SteadyJacobianPolicy & operator=(SteadyJacobianPolicy &&) = delete;
  ~SteadyJacobianPolicy() = default;

  SteadyJacobianPolicy(const FomSystemType & fomSystem,
		       FomStatesManagerType & fomStatesMngr,
		       DecoderType & decoder)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      decoderObj_(decoder){}

public:
  jacobian_type create() const
  {
    const auto & decJac = decoderObj_.get().jacobianCRef();
    return jacobian_type(fomSystem_.get().createApplyJacobianResult(decJac));
  }

  template <class lspg_state_t, class lspg_jac_t>
  void compute(const lspg_state_t & romState,
	       lspg_jac_t & romJacobian) const
  {
    // update Jacobian of decoder
    decoderObj_.get().updateJacobian(romState);

    // // todo: this is not needed if jacobian is called after resiudal
    // // because residual takes care of reconstructing the fom state
    // //    timer->start("reconstruct fom state");
    // fomStatesMngr_.get().template reconstructCurrentFomState(romState);

    const auto & decJac = decoderObj_.get().jacobianCRef();
    const auto & currFom = fomStatesMngr_.get().currentFomState();
    fomSystem_.get().applyJacobian(currFom, decJac, romJacobian);
  }

protected:
  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<DecoderType> decoderObj_;
};

}}}}
#endif  // ROM_IMPL_ROM_LSPG_STEADY_POLICIES_HPP_
