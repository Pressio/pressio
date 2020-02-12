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

#include "../rom_ConfigDefs.hpp"

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

    constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    // yOut = yOut + yFomReference_;
    containers::ops::do_update(yOut, one, yFomReference_, one);
  }

  template <typename rom_state_t>
  fom_state_type operator()(const rom_state_t & romY) const{
    auto yOut(yFomReference_);
    ::pressio::containers::ops::set_zero(yOut);
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  const fom_state_type & yFomReference_	= {};
  const decoder_type   & decoderObj_	= {};

};//end class


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <
  typename scalar_type,
  typename fom_state_type,
  typename decoder_type
  >
class FomStateReconstructorPressioOps<
  scalar_type, fom_state_type, decoder_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<fom_state_type>::value
    >
  >
{
public:
  FomStateReconstructorPressioOps() = delete;
  FomStateReconstructorPressioOps(const fom_state_type & yFomIn,
				  const decoder_type & decoder)
    : yFomReference_(yFomIn),
      decoderObj_(decoder)
  {}

public:
  template <typename rom_state_t>
  void operator()(const rom_state_t & romY,
		  fom_state_type    & yOut) const
  {
    decoderObj_.applyMapping(romY, yOut);
    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
    // add reference state yOut += yFomReference
    ::pressio::containers::ops::do_update(yOut, one, yFomReference_, one);
  }

  template <typename rom_state_t>
  fom_state_type operator()(const rom_state_t & romY) const{
    fom_state_type yOut{ fom_state_t(yFomReference_.request()) };
    ::pressio::containers::ops::set_zero(yOut);
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  fom_state_type yFomReference_	   = {};
  const decoder_type & decoderObj_ = {};

};//end class
#endif

}}}//end namespace pressio::rom::impl
#endif