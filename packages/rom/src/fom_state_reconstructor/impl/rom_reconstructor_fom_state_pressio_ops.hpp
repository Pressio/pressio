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

#ifndef ROM_RECONSTRUCTOR_FOM_STATE_PRESSIO_OPS_HPP_
#define ROM_RECONSTRUCTOR_FOM_STATE_PRESSIO_OPS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  typename scalar_type,
  typename fom_state_type,
  typename decoder_type
  >
struct FomStateReconstructorPressioOps
{
  FomStateReconstructorPressioOps() = delete;

  FomStateReconstructorPressioOps(const fom_state_type & yFomIn,
				  const decoder_type & decoder)
    : yFomReference_(yFomIn),
      decoderObj_(decoder)
  {}

  template <typename rom_state_t>
  void operator()(const rom_state_t & romY,
		  fom_state_type    & yOut) const
  {
    // map current romY to FOM state
    decoderObj_.applyMapping(romY, yOut);
    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    // yOut = yOut + yFomReference_;
    ops::do_update(yOut, one, yFomReference_, one);
  }

  template <typename rom_state_t>
  fom_state_type operator()(const rom_state_t & romY) const{
    auto yOut(yFomReference_);
    ::pressio::ops::set_zero(yOut);
    this->operator()(romY,yOut);
    return yOut;
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // this is added for pybind because I cannot figure out how to overload ()
  template <typename rom_state_t>
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_array_pybind<rom_state_t>::value,
    typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t
    >
  evaluate(const rom_state_t & romY) const{
    ::pressio::containers::Vector<rom_state_t> romView(romY, ::pressio::view());
    fom_state_type yOut(yFomReference_);
    ::pressio::ops::set_zero(yOut);
    this->operator()(romView, yOut);
    return *yOut.data();
  }
#endif

private:
  const fom_state_type & yFomReference_	= {};
  const decoder_type   & decoderObj_	= {};

};//end class

}}}//end namespace pressio::rom::impl
#endif
