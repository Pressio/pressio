/*
//@HEADER
// ************************************************************************
//
// rom_will_be_concepts.hpp
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

#ifndef ROM_WILL_BE_CONCEPTS_ROM_WILL_BE_CONCEPTS_INCLUDE_HPP_
#define ROM_WILL_BE_CONCEPTS_ROM_WILL_BE_CONCEPTS_INCLUDE_HPP_

// fom and rom states
#include "./various/rom_fom_state.hpp"
#include "./various/rom_rom_state.hpp"

// decoder
#include "./decoder/rom_admissible_decoder.hpp"
#include "./decoder/rom_decoder_jacobian.hpp"

// custom ops
#include "./custom_ops/rom_custom_ops_for_linear_decoder.hpp"
#include "./custom_ops/rom_custom_ops_for_fom_state_reconstructor.hpp"
#include "./custom_ops/rom_custom_ops_galerkin_continuous_time.hpp"
#include "./custom_ops/rom_custom_ops_lspg_continuous_time.hpp"
#include "./custom_ops/rom_custom_ops_lspg_discrete_time.hpp"

// masker and preconditioner
#include "./various/rom_masker.hpp"
#include "./various/rom_preconditioner.hpp"

// fom system
#include "./system/rom_steady_system.hpp"
#include "./system/rom_discrete_time_system_with_user_provided_apply_jacobian.hpp"
#include "./system/rom_continuous_time_system_without_user_provided_apply_jacobian.hpp"
#include "./system/rom_continuous_time_system_with_user_provided_apply_jacobian.hpp"
#include "./system/rom_continuous_time_system.hpp"

#endif  // ROM_WILL_BE_CONCEPTS_ROM_WILL_BE_CONCEPTS_INCLUDE_HPP_
