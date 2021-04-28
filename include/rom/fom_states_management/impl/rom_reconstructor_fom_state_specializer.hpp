/*
//@HEADER
// ************************************************************************
//
// rom_reconstructor_fom_state_specializer.hpp
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

#ifndef ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_SPECIALIZER_HPP_
#define ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_SPECIALIZER_HPP_

#include "rom_reconstructor_fom_state_pressio_ops.hpp"
#include "rom_reconstructor_fom_state_custom_ops.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <typename ... Args>
struct FomStateReconstructorSpecializer;

template <
  typename scalar_type,
  typename fom_state_type,
  typename decoder_type
  >
struct FomStateReconstructorSpecializer<
  scalar_type, fom_state_type, decoder_type
  >
{
  // check for a valid scalar type
  static_assert
  (std::is_floating_point<scalar_type>::value,
   "The first template arg to FomStateReconstructor must be a floating point type");

  // check for a valid FOM state type
  static_assert
  (::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
   "The second template arg to FomStateReconstructor must be a pressio wrapper.");

  using type = FomStateReconstructorPressioOps<scalar_type, fom_state_type, decoder_type>;
};

template <
  typename scalar_type,
  typename fom_state_type,
  typename decoder_type,
  typename ops_t
  >
struct FomStateReconstructorSpecializer<
  scalar_type, fom_state_type, decoder_type, ops_t
  >
{
  // check for a valid scalar type
  static_assert
  (std::is_floating_point<scalar_type>::value,
   "The first template arg to FomStateReconstructor must be a floating point type");

  // check for a valid FOM state type
  static_assert
  (::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
   "The second template arg to FomStateReconstructor must be a pressio wrapper.");

  // check for valid ops
  static_assert
  (::pressio::rom::constraints::custom_ops_for_fom_state_reconstructor<
   ops_t, fom_state_type >::value,
   "You are tring to create a FomStateReconstructor with custom ops. \
The template arg passed representing the custom ops \
does not have an admissible API for the FomStateReconstructor");

  using type = FomStateReconstructorCustomOps<scalar_type, fom_state_type, decoder_type, ops_t>;
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_FOM_STATES_MANAGEMENT_IMPL_ROM_RECONSTRUCTOR_FOM_STATE_SPECIALIZER_HPP_
