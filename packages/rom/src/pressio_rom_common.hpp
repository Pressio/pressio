/*
//@HEADER
// ************************************************************************
//
// pressio_rom.hpp
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

#ifndef PRESSIO_ROM_COMMON_HPP_
#define PRESSIO_ROM_COMMON_HPP_

/*
   NOTE that the order below matters!
   The clean/logical order allows us to avoid ending up with a tangled system.
   NOTE also that this header by itself means nothing and if you use
   it as such, you need to know what you are doing.
   This header is placed here to help the "public" includes
   named "pressio_rom_bla.hpp" at the top level.
   Users of pressio should NOT rely on this, that is why this is
   placed here. Users should rely only on the top-level
   "pressio_rom_{lspg,galerkin,wls}.hpp".
*/

// need forward declarations
#include "rom_fwd.hpp"

// all predicates
#include "./predicates/rom_predicates.hpp"

// all will_be_concepts (depend on predicates)
#include "./will_be_concepts/rom_will_be_concepts.hpp"

// decoder classes (depend on concepts)
#include "./decoder/rom_decoders.hpp"

// fom states management classes (depend on the decoder)
#include "./fom_states_management/rom_manager_fom_states_static.hpp"
#include "./fom_states_management/rom_reconstructor_fom_state.hpp"

// decorators
#include "./decorators/rom_decorators.hpp"

#endif
