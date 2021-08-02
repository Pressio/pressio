/*
//@HEADER
// ************************************************************************
//
// pressio_solvers.hpp
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

#ifndef PRESSIO_NONLINEAR_SOLVERS_HPP_
#define PRESSIO_NONLINEAR_SOLVERS_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_type_traits.hpp"
#include "pressio_expressions.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"

#include "solvers_nonlinear/solvers_exceptions.hpp"

#include "solvers_nonlinear/solvers_nonlinear_tags.hpp"
#include "solvers_nonlinear/solvers_nonlinear_enums.hpp"
#include "solvers_nonlinear/constraints/solvers_implicit_state.hpp"
#include "solvers_nonlinear/constraints/solvers_least_squares_weighting_operator.hpp"
#include "solvers_nonlinear/impl/updaters/solvers_create_updater.hpp"
#include "solvers_nonlinear/impl/solvers_observer.hpp"

#include "solvers_nonlinear/constraints/solvers_legitimate_linear_solver_for_newton_raphson.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_residualnorm_method_accept_state_norm_return_void.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_create_residual_method_return_result.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_create_jacobian_method_return_result.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_jacobian_method_accept_state_result_return_void.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_residual_method_accept_state_result_return_void.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_residualandjacobian_method_accept_state_result_return_void.hpp"
#include "solvers_nonlinear/constraints/solvers_system_fused_residual_jacobian.hpp"
#include "solvers_nonlinear/constraints/solvers_system_residual_jacobian.hpp"
#include "solvers_nonlinear/impl/solvers_printer.hpp"

// *** non-linear least-squares *** //
#include "solvers_nonlinear/predicates/solvers_has_const_create_hessian_method_return_result.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_create_gradient_method_return_result.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_gradient_method_accept_state_result_norm_return_void.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_hessian_method_accept_state_result_return_void.hpp"
#include "solvers_nonlinear/predicates/solvers_has_const_hessianandgradient_method_accept_state_result_norm_return_void.hpp"
#include "solvers_nonlinear/constraints/solvers_system_fused_hessian_gradient.hpp"
#include "solvers_nonlinear/constraints/solvers_system_hessian_gradient.hpp"
#include "solvers_nonlinear/constraints/solvers_legitimate_linear_solver_for_nonlinear_least_squares.hpp"
#include "solvers_nonlinear/constraints/solvers_legitimate_qr_solver_for_gn_qr.hpp"

#include "solvers_nonlinear/solvers_create_public_api.hpp"

#endif