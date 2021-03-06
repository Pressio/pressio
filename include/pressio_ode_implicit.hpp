/*
//@HEADER
// ************************************************************************
//
// pressio_ode_implicit.hpp
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

#ifndef PRESSIO_ODE_IMPLICIT_HPP_
#define PRESSIO_ODE_IMPLICIT_HPP_

/*
   include everything needed for ODE implicit integration
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
#include "pressio_solvers.hpp"

// common things
#include "ode/pressio_ode_common.hpp"

// specific to implicit
#include "ode/implicit/ode_implicit_constants.hpp"
#include "ode/implicit/constraints/ode_implicit_state.hpp"
#include "ode/implicit/constraints/ode_implicit_residual.hpp"
#include "ode/implicit/constraints/ode_implicit_jacobian.hpp"
#include "ode/implicit/constraints/ode_legitimate_solver_for_implicit_stepper.hpp"
#include "ode/implicit/constraints/ode_user_defined_ops_for_implicit_bdf2.hpp"
#include "ode/implicit/constraints/ode_user_defined_ops_for_implicit_euler.hpp"
#include "ode/implicit/constraints/ode_user_defined_ops_for_implicit_ode.hpp"
#include "ode/implicit/constraints/ode_implicitly_steppable.hpp"
#include "ode/implicit/constraints/ode_implicitly_steppable_with_guesser.hpp"
#include "ode/implicit/constraints/ode_auxiliary_stepper_for_bdf2.hpp"

#include "ode/implicit/ode_stencil_states_manager.hpp"
#include "ode/implicit/ode_stencil_velocities_manager.hpp"
#include "ode/implicit/constraints/ode_implicit_residual_policy.hpp"
#include "ode/implicit/constraints/ode_implicit_jacobian_policy.hpp"

#include "ode/implicit/ode_implicit_stepper.hpp"

#include "ode/integrators/ode_advance_n_steps_implicit_arbitrary_step_size.hpp"
#include "ode/integrators/ode_advance_n_steps_implicit_constant_step_size.hpp"
#include "ode/integrators/ode_advance_to_target_time_implicit_arbitrary_step_size.hpp"

#endif
