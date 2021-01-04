/*
//@HEADER
// ************************************************************************
//
// pressio_ode_common.hpp
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

#ifndef ODE_PRESSIO_ODE_COMMON_HPP_
#define ODE_PRESSIO_ODE_COMMON_HPP_

/*
   NOTE that the order below matters!
   - Includes are ordered properly to avoid a tangled system.
   - don't rely on files inside impl, these might change

   NOTE also that this header by itself means nothing and if you use
   it as such, you need to know what you are doing.
   This header is here to help the publicly-exposed includes named
   "pressio_ode_bla.hpp" inside the pressio/packages directory.
   Users of pressio should NOT rely on this file, but only
   on the top-level "pressio_ode_{explicit,implicit}.hpp".
*/

//----------------------------------------------
#include "ode/src/ode_types.hpp"
#include "ode/src/ode_stencil_tags.hpp"
#include "ode/src/ode_exceptions.hpp"
#include "ode/src/ode_stepper_tags.hpp"
#include "ode/src/ode_required_number_of_states.hpp"

//----------------------------------------------
// predicates
#include "ode/src/predicates/typedefs/ode_has_state_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_velocity_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_residual_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_jacobian_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_discrete_time_residual_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_discrete_time_jacobian_typedef.hpp"

#include "ode/src/predicates/velocity_methods/ode_has_const_create_velocity_method_return_result.hpp"
#include "ode/src/predicates/velocity_methods/ode_has_const_velocity_method_accept_state_time_result_return_void.hpp"
#include "ode/src/predicates/ode_is_stepper_total_n_states_setter.hpp"
#include "ode/src/predicates/ode_is_stepper_order_setter.hpp"
#include "ode/src/predicates/discrete_time_residual_methods/ode_has_const_create_discrete_time_residual_method_return_result.hpp"
#include "ode/src/predicates/discrete_time_residual_methods/ode_has_const_discrete_time_residual_method_accept_step_time_dt_result_states_return_void.hpp"
#include "ode/src/predicates/discrete_time_jacobian_methods/ode_has_const_create_discrete_time_jacobian_method_return_result.hpp"
#include "ode/src/predicates/discrete_time_jacobian_methods/ode_has_const_discrete_time_jacobian_method_accepting_n_states_returning_void.hpp"
#include "ode/src/predicates/jacobian_methods/ode_has_const_create_jacobian_method_return_result.hpp"
#include "ode/src/predicates/jacobian_methods/ode_has_const_jacobian_method_accept_state_time_result_return_void.hpp"

//----------------------------------------------
// common constraints (depend on predicates)
#include "ode/src/constraints/ode_collector.hpp"
#include "ode/src/constraints/ode_guesser.hpp"
#include "ode/src/constraints/ode_time_step_size_manager.hpp"
#include "ode/src/constraints/system/ode_continuous_time_system_with_at_least_velocity.hpp"
#include "ode/src/constraints/system/ode_continuous_time_system_with_user_provided_jacobian.hpp"
#include "ode/src/constraints/system/ode_discrete_time_system_with_user_provided_jacobian.hpp"

#endif  // ODE_PRESSIO_ODE_COMMON_HPP_
