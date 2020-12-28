/*
//@HEADER
// ************************************************************************
//
// pressio_rom_galerkin.hpp
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

#ifndef PRESSIO_ROM_GALERKIN_HPP_
#define PRESSIO_ROM_GALERKIN_HPP_

/*
   this header includes everything needed for GALERKIN.
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
#include "pressio_ode.hpp"

// common classes for rom
#include "rom/src/pressio_rom_common.hpp"

// constraints
#include "rom/src/galerkin/constraints/rom_fom_state.hpp"
#include "rom/src/galerkin/constraints/rom_fom_velocity.hpp"

#include "rom/src/galerkin/constraints/rom_galerkin_explicit_state.hpp"
#include "rom/src/galerkin/constraints/rom_galerkin_implicit_state.hpp"
#include "rom/src/galerkin/constraints/rom_galerkin_velocity.hpp"
#include "rom/src/galerkin/constraints/rom_galerkin_residual.hpp"
#include "rom/src/galerkin/constraints/rom_galerkin_jacobian.hpp"

#include "rom/src/galerkin/constraints/rom_custom_ops_continuous_time.hpp"
#include "rom/src/galerkin/constraints/rom_masker_explicit.hpp"
#include "rom/src/galerkin/constraints/rom_masker_implicit.hpp"
#include "rom/src/galerkin/constraints/rom_projector_explicit.hpp"
#include "rom/src/galerkin/constraints/rom_projector_implicit.hpp"

// projectors
#include "rom/src/galerkin/impl/projectors/galerkin_arbitrary_projector.hpp"
#include "rom/src/galerkin/impl/projectors/galerkin_default_projector_ortho_decoder_jacobian.hpp"

// galerkin classes
#include "rom/src/impl/rom_auxiliary_stepper_type_helper.hpp"
#include "rom/src/impl/rom_problem_members_mixins.hpp"
#include "rom/src/galerkin/impl/rom_problem_members.hpp"

#include "rom/src/galerkin/impl/decorators/rom_masked.hpp"
#include "rom/src/galerkin/impl/decorators/rom_projected.hpp"
#include "rom/src/galerkin/impl/policies/rom_galerkin_residual_policy.hpp"
#include "rom/src/galerkin/impl/policies/rom_galerkin_jacobian_policy.hpp"

#include "rom/src/galerkin/impl/continuous_time_api/rom_compose_impl.hpp"
#include "rom/src/galerkin/impl/discrete_time_api/rom_compose_impl.hpp"

#include "rom/src/galerkin/rom_create_default_galerkin_problem.hpp"
#include "rom/src/galerkin/rom_create_masked_galerkin_problem.hpp"
#include "rom/src/galerkin/rom_create_hyperreduced_galerkin_problem.hpp"

#include "rom/src/galerkin/rom_create_solver_functions.hpp"
#include "rom/src/galerkin/rom_solve_problem_functions.hpp"
#include "rom/src/galerkin/rom_create_collocation_projector.hpp"
#include "rom/src/galerkin/rom_create_arbitrary_projector.hpp"

#endif
