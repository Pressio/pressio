/*
//@HEADER
// ************************************************************************
//
// pressio_rom_lspg.hpp
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

#ifndef PRESSIO_ROM_LSPG_INCLUDE_HPP_
#define PRESSIO_ROM_LSPG_INCLUDE_HPP_

/*
   this header includes everything needed for LSPG.
   NOTE that the order below matters!
   - Includes are ordered properly to avoid a tangled system.
   - don't rely on files inside impl, these might change
*/

// need all of the dependent packages
#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"
#include "pressio_optimizers.hpp"
#include "pressio_solvers.hpp"
#include "pressio_ode_implicit.hpp"

// common classes for rom
#include "rom/src/pressio_rom_common.hpp"

// constraints
#include "rom/src/lspg/constraints/rom_fom_state.hpp"
#include "rom/src/lspg/constraints/rom_fom_velocity.hpp"
#include "rom/src/lspg/constraints/rom_lspg_jacobian.hpp"
#include "rom/src/lspg/constraints/rom_lspg_residual.hpp"
#include "rom/src/lspg/constraints/rom_lspg_state.hpp"
#include "rom/src/lspg/constraints/rom_steady_preconditioner.hpp"
#include "rom/src/lspg/constraints/rom_unsteady_preconditioner.hpp"
#include "rom/src/lspg/constraints/rom_steady_masker.hpp"
#include "rom/src/lspg/constraints/rom_unsteady_masker.hpp"
#include "rom/src/lspg/constraints/rom_custom_ops_continuous_time.hpp"
#include "rom/src/lspg/constraints/rom_custom_ops_discrete_time.hpp"

// decorators
#include "rom/src/lspg/decorators/rom_preconditioned.hpp"
#include "rom/src/lspg/decorators/rom_masked.hpp"

// lspg classes
#include "rom/src/impl/rom_auxiliary_stepper_type_helper.hpp"
#include "rom/src/impl/rom_problem_members_mixins.hpp"
#include "rom/src/lspg/impl/rom_problem_members.hpp"
#include "rom/src/lspg/impl/steady/rom_compose_steady_lspg_impl.hpp"
#include "rom/src/lspg/impl/unsteady/rom_compose_unsteady_lspg_impl.hpp"
#include "rom/src/lspg/impl/rom_compose_problem.hpp"

#include "rom/src/lspg/rom_create_default_lspg_problem.hpp"
#include "rom/src/lspg/rom_create_preconditioned_default_lspg_problem.hpp"
#include "rom/src/lspg/rom_create_masked_lspg_problem.hpp"
#include "rom/src/lspg/rom_create_hyper_reduced_lspg_problem.hpp"
#include "rom/src/lspg/rom_create_preconditioned_hyper_reduced_lspg_problem.hpp"
#include "rom/src/lspg/rom_create_solver_functions.hpp"
#include "rom/src/lspg/rom_solve_problem_functions.hpp"

#endif
