/*
//@HEADER
// ************************************************************************
//
// rom_reconstructor_fom_state.hpp
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

#ifndef ROM_RECONSTRUCTOR_FOM_STATE_HPP_
#define ROM_RECONSTRUCTOR_FOM_STATE_HPP_

#include "rom_ConfigDefs.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_state_t,
  typename decoder_type,
  typename enable = void
  >
struct FomStateReconstructor;


template <
  typename fom_state_t,
  typename decoder_type
  >
struct FomStateReconstructor<
  fom_state_t, decoder_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<fom_state_t>::value
    >
  >
{
  using scalar_t = typename containers::details::traits<fom_state_t>::scalar_t;

  FomStateReconstructor() = delete;

  FomStateReconstructor(const fom_state_t & yFomIn,
			const decoder_type & decoder)
    : yFomReference_(yFomIn),
      decoderObj_(decoder)
  {}

  ~FomStateReconstructor() = default;

  template <typename rom_state_t>
  void operator()(const rom_state_t	& romY,
		  fom_state_t		& yOut) const
  {
    // map current romY to FOM state
    decoderObj_.applyMapping(romY, yOut);

    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
    // yOut = yOut + yFomReference_;
    containers::ops::do_update(yOut, one, yFomReference_, one);
  }

  template <typename rom_state_t>
  fom_state_t operator()(const rom_state_t & romY) const{
    auto yOut(yFomReference_);
    ::pressio::containers::ops::set_zero(yOut);
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  const fom_state_t & yFomReference_	= {};
  const decoder_type & decoderObj_	= {};

};//end class






#ifdef HAVE_PYBIND11
template <
  typename fom_state_t,
  typename decoder_type
  >
class FomStateReconstructor<
  fom_state_t, decoder_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<fom_state_t>::value
    >
  >
{
  using scalar_t = typename fom_state_t::value_type;

public:
  FomStateReconstructor() = delete;
  ~FomStateReconstructor() = default;

  FomStateReconstructor(const fom_state_t & yFomIn,
			const decoder_type & decoder)
    : /*yFomReference_( const_cast<fom_state_t&>(yFomIn).request().shape,
		      const_cast<fom_state_t&>(yFomIn).request().strides,
		      yFomIn.data() ),*/
      yFomReference_(yFomIn),
      decoderObj_(decoder)
  {
#ifdef DEBUG_PRINT
    std::cout << std::endl;
    std::cout << "Inside FomStateReconstructor " << std::endl;
    std::cout << "yFomReference_ " << yFomReference_.data() << std::endl;
    std::cout << std::endl;
#endif
  }

public:
  template <typename rom_state_t>
  void operator()(const rom_state_t	& romY,
		  fom_state_t		& yOut) const {
    decoderObj_.applyMapping(romY, yOut);

    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
    ::pressio::containers::ops::do_update(yOut, one, yFomReference_, one);
    //yOut += yFomReference_;
  }

  template <typename rom_state_t>
  fom_state_t operator()(const rom_state_t & romY) const{
    fom_state_t yOut{ fom_state_t(yFomReference_.request()) };
    ::pressio::containers::ops::set_zero(yOut);
    this->template operator()(romY,yOut);
    return yOut;
  }

private:
  fom_state_t yFomReference_	= {};
  const decoder_type & decoderObj_	= {};

};//end class
#endif

}}//end namespace pressio::rom
#endif
