/*
//@HEADER
// ************************************************************************
//
// rom_reconstructor_fom_state_pressio_ops.hpp
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

#ifndef ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_HPP_
#define ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class DecoderType>
struct FomStateReconstructor
{
  using decoder_type = mpl::remove_cvref_t<DecoderType>;

  using fom_state_type = typename decoder_type::fom_state_type;
  using scalar_type = typename ::pressio::Traits<fom_state_type>::scalar_type;

  FomStateReconstructor() = delete;
  FomStateReconstructor(const FomStateReconstructor &) = default;
  FomStateReconstructor & operator=(const FomStateReconstructor &) = delete;
  FomStateReconstructor(FomStateReconstructor &&) = default;
  FomStateReconstructor & operator=(FomStateReconstructor &&) = delete;
  ~FomStateReconstructor() = default;

  FomStateReconstructor(const fom_state_type & fomNominalState,
			const decoder_type & decoder)
    : fomNominalState_(::pressio::ops::clone(fomNominalState)),
      decoderObj_(decoder)
  {}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   template <typename rom_state_t>
//   void operator()(const rom_state_t & romState,
// 		  fom_state_type    & fomState) const
//   {
//     // map current romState to FOM state
//     decoderObj_.get().applyMapping(romState, fomState);
//     constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
//     // fomState = fomState + fomNominalState_;
//     ops::update(fomState, one, fomNominalState_, one);
//   }

//   // evaluate is added because I cannot figure out how to overload ()
//    evaluate is used to pass an array from python and construct
//      the fom result here which we pass back to pyhon.
//      So the rom state is a native python array.

//      Note that we template this on the rom_state_t but we pass the native array.
//      By doing this, the template allows us to know the rank of the rom state
//      since this should work for rank1 and rank2 rom states.
//      If we only passed the native array wihtout the template, we would not know
//      the propper rank of the fom state.

//   template<typename rom_state_t>
//   native_fom_state_t evaluate(const typename rom_state_t::traits::wrapped_t & romStateIn) const
//   {
//     rom_state_t romView(romStateIn, ::pressio::view());
//     fom_state_type fomState(fomNominalState_);
//     ::pressio::ops::set_zero(fomState);
//     (*this)(romView, fomState);
//     return *fomState.data();
//   }

// #else
  template <class RomStateType>
  void operator()(const RomStateType & romState,
		  fom_state_type & fomState) const
  {

    // map current romState to FOM state
    decoderObj_.get().applyMapping(romState, fomState);
    constexpr auto one = ::pressio::utils::Constants<scalar_type>::one();
    // fomState = fomState + fomNominalState_;
    ::pressio::ops::update(fomState, one, fomNominalState_, one);
  }

  template <class RomStateType>
  fom_state_type operator()(const RomStateType & romState) const
  {

    auto fomState = ::pressio::ops::clone(fomNominalState_);
    ::pressio::ops::set_zero(fomState);
    this->operator()(romState,fomState);
    return fomState;
  }
// #endif

private:
  const fom_state_type fomNominalState_	= {};
  std::reference_wrapper<const decoder_type> decoderObj_ = {};
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_PRESSIO_OPS_HPP_