/*
//@HEADER
// ************************************************************************
//
// rom_wls_select_timescheme_helper.hpp
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

#ifndef ROM_WLS_TIME_SCHEMES_ROM_WLS_SELECT_TIMESCHEME_HELPER_HPP_
#define ROM_WLS_TIME_SCHEMES_ROM_WLS_SELECT_TIMESCHEME_HELPER_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{

namespace impl{

template <typename ode_tag, typename fom_state_t>
struct DefaultHelper{
  using type = void;
};

//BDF1 specialization
template <typename fom_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::Euler, fom_state_t>{
  using type = ::pressio::rom::wls::timeschemes::impl::ImplicitEuler<fom_state_t>;
};

//BDF2 specialization
template <typename fom_state_t>
struct DefaultHelper<::pressio::ode::implicitmethods::BDF2, fom_state_t>{
  using type = ::pressio::rom::wls::timeschemes::impl::BDF2<fom_state_t>;
};

} //end namespace pressio::rom::wls::timeschemes::helpers::impl

template <typename ode_tag, typename fom_state_t>
using timescheme_t = typename impl::DefaultHelper<ode_tag,fom_state_t>::type;

}}}} //ena namespace pressio::rom::wls::timeschemes::helpers
#endif  // ROM_WLS_TIME_SCHEMES_ROM_WLS_SELECT_TIMESCHEME_HELPER_HPP_
