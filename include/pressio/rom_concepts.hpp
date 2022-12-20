/*
//@HEADER
// ************************************************************************
//
// rom_galerkin.hpp
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

#ifndef PRESSIO_ROM_CONCEPTS_TOPLEVEL_INC_HPP_
#define PRESSIO_ROM_CONCEPTS_TOPLEVEL_INC_HPP_

#include "./mpl.hpp"
#include "./utils.hpp"
#include "./type_traits.hpp"
#include "./ops.hpp"
#include "./qr.hpp"
#include "./solvers_linear.hpp"
#include "./solvers_nonlinear.hpp"
#include "./ode.hpp"

#ifdef PRESSIO_ENABLE_CXX20
#include <concepts>
#endif

#include "./rom/predicates.hpp"
#include "./rom/reduced_operators_traits.hpp"
#include "./rom/concepts/valid_reduced_state.hpp"

#include "./rom/concepts/maskable_with.hpp"

// subspaces
#include "./rom/concepts/linear_subspace.hpp"
#include "./rom/concepts/possibly_affine_trial_column_subspace.hpp"

// fom
#include "./rom/concepts/fom_steady_with_jac_action.hpp"
#include "./rom/concepts/fom_semi_discrete.hpp"
#include "./rom/concepts/fom_semi_discrete_with_mm_action.hpp"
#include "./rom/concepts/fom_semi_discrete_with_jac_action.hpp"
#include "./rom/concepts/fom_fully_discrete_with_jac_action.hpp"
#include "./rom/concepts/fom_semi_discrete_with_jac_and_mm_action.hpp"

// galerkin
#include "./rom/concepts/galerkin_steady_default.hpp"
#include "./rom/concepts/galerkin_steady_hyperreduced.hpp"
#include "./rom/concepts/galerkin_steady_masked.hpp"
#include "./rom/concepts/galerkin_explicit_default.hpp"
#include "./rom/concepts/galerkin_explicit_default_with_varying_mm.hpp"
#include "./rom/concepts/galerkin_explicit_hyperreduced.hpp"
#include "./rom/concepts/galerkin_explicit_masked.hpp"
#include "./rom/concepts/galerkin_implicit_default.hpp"
#include "./rom/concepts/galerkin_implicit_default_with_varying_mm.hpp"
#include "./rom/concepts/galerkin_implicit_hyperreduced.hpp"
#include "./rom/concepts/galerkin_implicit_masked.hpp"


// lspg
#include "./rom/concepts/lspg_steady_default.hpp"
#include "./rom/concepts/lspg_steady_masked.hpp"
#include "./rom/concepts/lspg_unsteady_default.hpp"
#include "./rom/concepts/lspg_unsteady_hyperreduced.hpp"
#include "./rom/concepts/lspg_unsteady_masked.hpp"

#endif
