/*
//@HEADER
// ************************************************************************
//
// rom_reconstructor_fom_state_custom_ops.hpp
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

#ifndef ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_CUSTOM_OPS_HPP_
#define ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_CUSTOM_OPS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  typename scalar_type,
  typename fom_state_type,
  typename decoder_type,
  typename ops_type
  >
struct FomStateReconstructorCustomOps
{
  FomStateReconstructorCustomOps() = delete;
  FomStateReconstructorCustomOps(const FomStateReconstructorCustomOps &) = default;
  FomStateReconstructorCustomOps & operator=(const FomStateReconstructorCustomOps &) = delete;
  FomStateReconstructorCustomOps(FomStateReconstructorCustomOps &&) = default;
  FomStateReconstructorCustomOps & operator=(FomStateReconstructorCustomOps &&) = delete;
  ~FomStateReconstructorCustomOps() = default;

  FomStateReconstructorCustomOps(const fom_state_type & fomStateIn,
				 const decoder_type & decoder,
				 const ops_type & udOps)
    : fomNominalState_(fomStateIn), decoderObj_(decoder), udOps_{udOps}
  {
    PRESSIOLOG_DEBUG("cnstr: fomNominalState extent = {}",
		     fomStateIn.extent(0));
  }

  template <typename rom_state_t>
  void operator()(const rom_state_t & romState,
		  fom_state_type    & fomState) const
  {
    // map current romState to FOM state
    decoderObj_.get().applyMapping(romState, fomState);

    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    // fomState = fomState + fomNominalState_;
    udOps_.get().axpy(one, *(fomNominalState_.get().data()), *fomState.data());
  }

  template <typename rom_state_t>
  fom_state_type operator()(const rom_state_t & romState) const{
    auto fomState(fomNominalState_.get());
    udOps_.get().set_zero(*fomState.data());
    this->template operator()(romState,fomState);
    return fomState;
  }

private:
  std::reference_wrapper<const fom_state_type> fomNominalState_;
  std::reference_wrapper<const decoder_type> decoderObj_;
  std::reference_wrapper<const ops_type> udOps_;

};//end class

}}}//end namespace pressio::rom::impl
#endif  // ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_CUSTOM_OPS_HPP_
