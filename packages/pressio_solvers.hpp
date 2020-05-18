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

#ifndef PRESSIO_SOLVERS_HPP_
#define PRESSIO_SOLVERS_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"

#include "solvers/src/solvers_norm_tags.hpp"
#include "solvers/src/solvers_convergence_tags.hpp"
#include "solvers/src/solvers_fwd.hpp"
#include "solvers/src/solvers_meta_static_checks.hpp"

// base classes
#include "solvers/src/base/solvers_iterative_base.hpp"
#include "solvers/src/base/solvers_linear_base.hpp"
#include "solvers/src/base/solvers_nonlinear_base.hpp"

// meta
#include "solvers/src/meta/solvers_basic_meta.hpp"


//**********************
// *** linear *** //
//**********************
#include "solvers/src/linear/solvers_linear_tags.hpp"
#include "solvers/src/linear/solvers_linear_traits.hpp"
#include "solvers/src/linear/solvers_linear_solver.hpp"


//**********************
// *** non-linear *** //
//**********************
#include "solvers/src/meta/custom_ops_detection/solvers_has_all_needed_methods_norms.hpp"
#include "solvers/src/meta/custom_ops_detection/solvers_has_all_needed_methods_for_hessian.hpp"
#include "solvers/src/meta/custom_ops_detection/solvers_has_all_needed_methods_for_gradient.hpp"

#include "solvers/src/nonlinear/solvers_line_search_tags.hpp"
#include "solvers/src/meta/solvers_is_legitimate_convergence_tag.hpp"
#include "solvers/src/meta/solvers_is_legitimate_line_search_tag.hpp"
#include "solvers/src/meta/solvers_is_legitimate_hessian_for_gn_normeq.hpp"
#include "solvers/src/meta/solvers_is_legitimate_linear_solver_for_gn_normeq.hpp"
#include "solvers/src/meta/solvers_is_legitimate_qr_solver_for_gn_qr.hpp"
#include "solvers/src/meta/solvers_is_legitimate_residual_observer_each_solver_step.hpp"
#include "solvers/src/meta/solvers_is_legitimate_residual_observer_when_converged.hpp"
#include "solvers/src/meta/solvers_system_meets_gn_hessian_gradient_api.hpp"

#include "solvers/src/meta/solvers_system_has_all_needed_jacobian_methods.hpp"
#include "solvers/src/meta/solvers_system_has_all_needed_residual_methods.hpp"
#include "solvers/src/meta/solvers_system_meets_default_api.hpp"
#include "solvers/src/meta/solvers_is_legitimate_system_for_gauss_newton.hpp"
#include "solvers/src/meta/solvers_is_legitimate_system_for_newton_raphson.hpp"

// newton-raphson
#include "solvers/src/nonlinear/solvers_newton_raphson.hpp"

// Gauss-Newton
#include "solvers/src/nonlinear/gauss_newton/solvers_gauss_newton.hpp"

// Gauss-Newton conservative rom
#include "solvers/src/nonlinear/gn_conservative_rom/solvers_gauss_newton_conservative.hpp"

#endif
