/*
//@HEADER
// ************************************************************************
//
// solvers_convergence_criteria.hpp
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

#ifndef PRESSIO_SOLVERS_CONVERGENCE_CRITERIA_HPP_
#define PRESSIO_SOLVERS_CONVERGENCE_CRITERIA_HPP_

#include "./impl/convergence_mixins/solvers_stop_at_max_iters.hpp"
#include "./impl/convergence_mixins/solvers_stop_when_correction_norm_below_tol.hpp"
#include "./impl/convergence_mixins/solvers_stop_when_residual_norm_below_tol.hpp"
#include "./impl/convergence_mixins/solvers_stop_when_gradient_norm_below_tol.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{

// correction norm below tol
template<typename ... Args>
using ConvergedWhenCorrectionNormBelowTol = impl::StopWhenCorrectionNormBelowTol<Args...>;
template<typename ... Args>
using StopWhenCorrectionNormBelowTol      = impl::StopWhenCorrectionNormBelowTol<Args...>;

// default convergence
template<typename ... Args>
using DefaultConvergence = impl::StopWhenCorrectionNormBelowTol<Args...>;

//--------------------------------
// stop at max iters
template<typename ... Args> using ConvergedWhenMaxIters = impl::StopAtMaxIters<Args...>;
template<typename ... Args> using StopWhenMaxIters	= impl::StopAtMaxIters<Args...>;
template<typename ... Args> using StopAfterMaxIters	= impl::StopAtMaxIters<Args...>;
template<typename ... Args> using StopAtMaxIters	= impl::StopAtMaxIters<Args...>;

//--------------------------------
// *** residual norm ***
//--------------------------------
// absolute norm below tol
template<typename ... Args>
using ConvergedWhenAbsoluteResidualNormBelowTol = impl::StopWhenResidualNormBelowTol<true, Args...>;
template<typename ... Args>
using StopWhenAbsoluteResidualNormBelowTol	= impl::StopWhenResidualNormBelowTol<true, Args...>;

// relative residual norm below tol
template<typename ... Args>
using ConvergedWhenRelativeResidualNormBelowTol = impl::StopWhenResidualNormBelowTol<false, Args...>;
template<typename ... Args>
using StopWhenRelativeResidualNormBelowTol	= impl::StopWhenResidualNormBelowTol<false, Args...>;

//--------------------------------
// *** gradiet norm ***
//--------------------------------
// absolute gradient norm below tol
template<typename ... Args>
using ConvergedWhenAbsoluteGradientNormBelowTol = impl::StopWhenGradientNormBelowTol<true, Args...>;
template<typename ... Args>
using StopWhenAbsoluteGradientNormBelowTol	= impl::StopWhenGradientNormBelowTol<true, Args...>;

// relative gradient norm below tol
template<typename ... Args>
using ConvergedWhenRelativeGradientNormBelowTol = impl::StopWhenGradientNormBelowTol<false, Args...>;
template<typename ... Args>
using StopWhenRelativeGradientNormBelowTol	= impl::StopWhenGradientNormBelowTol<false, Args...>;

}}}
#endif
