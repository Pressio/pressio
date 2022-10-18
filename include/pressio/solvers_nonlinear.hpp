/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear.hpp
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

#include "./mpl.hpp"
#include "./utils.hpp"
#include "./type_traits.hpp"
#include "./expressions.hpp"
#include "./ops.hpp"
#include "./qr.hpp"

#include "solvers_nonlinear/solvers_exceptions.hpp"

#include "solvers_nonlinear/solvers_nonlinear_enums_and_tags.hpp"

#include "solvers_nonlinear/concepts/solvers_predicates.hpp"
#include "solvers_nonlinear/concepts/solvers_system_residual_jacobian.hpp"
#include "solvers_nonlinear/concepts/solvers_system_fused_residual_jacobian.hpp"
#include "solvers_nonlinear/concepts/solvers_system_hessian_gradient.hpp"
#include "solvers_nonlinear/concepts/solvers_system_fused_hessian_gradient.hpp"
#include "solvers_nonlinear/concepts/solvers_least_squares_weighting_operator.hpp"
#include "solvers_nonlinear/concepts/solvers_linear_solver_for_newton_raphson.hpp"
#include "solvers_nonlinear/concepts/solvers_linear_solver_for_nonlinear_least_squares.hpp"
#include "solvers_nonlinear/concepts/solvers_qr_solver_for_gn_qr.hpp"

#include "solvers_nonlinear/impl/updaters/solvers_create_updater.hpp"
#include "solvers_nonlinear/impl/solvers_observer.hpp"
#include "solvers_nonlinear/impl/solvers_printer.hpp"
#include "solvers_nonlinear/solvers_create_public_api.hpp"

#endif
