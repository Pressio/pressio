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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CONVERGENCE_CRITERIA_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CONVERGENCE_CRITERIA_HPP_

#include "./impl/solve_until_mixins/solvers_solve_until_max_iters.hpp"
#include "./impl/solve_until_mixins/solvers_solve_until_correction_norm_below_tol.hpp"
#include "./impl/solve_until_mixins/solvers_solve_until_residual_norm_below_tol.hpp"
#include "./impl/solve_until_mixins/solvers_solve_until_gradient_norm_below_tol.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{

// default convergence
template<typename ... Args>
using DefaultConvergence = impl::SolveUntilCorrectionNormBelowTol<Args...>;

//-------------------------------
// correction norm below tol
//-------------------------------
template<typename ... Args>
using ConvergedWhenCorrectionNormBelowTol = impl::SolveUntilCorrectionNormBelowTol<Args...>;
template<typename ... Args>
using StopWhenCorrectionNormBelowTol      = impl::SolveUntilCorrectionNormBelowTol<Args...>;

//--------------------------------
// stop at max iters
//--------------------------------
template<typename ... Args> using ConvergedWhenMaxIters = impl::SolveUntilMaxIters<Args...>;
template<typename ... Args> using StopWhenMaxIters	= impl::SolveUntilMaxIters<Args...>;
template<typename ... Args> using StopAfterMaxIters	= impl::SolveUntilMaxIters<Args...>;
template<typename ... Args> using StopAtMaxIters	= impl::SolveUntilMaxIters<Args...>;

//--------------------------------
// residual norm
//--------------------------------
// absolute norm below tol
template<typename ... Args>
using ConvergedWhenAbsoluteResidualNormBelowTol = impl::SolveUntilResidualNormBelowTol<true, Args...>;
template<typename ... Args>
using StopWhenAbsoluteResidualNormBelowTol	= impl::SolveUntilResidualNormBelowTol<true, Args...>;

// relative residual norm below tol
template<typename ... Args>
using ConvergedWhenRelativeResidualNormBelowTol = impl::SolveUntilResidualNormBelowTol<false, Args...>;
template<typename ... Args>
using StopWhenRelativeResidualNormBelowTol	= impl::SolveUntilResidualNormBelowTol<false, Args...>;

//--------------------------------
// gradiet norm
//--------------------------------
// absolute gradient norm below tol
template<typename ... Args>
using ConvergedWhenAbsoluteGradientNormBelowTol = impl::SolveUntilGradientNormBelowTol<true, Args...>;
template<typename ... Args>
using StopWhenAbsoluteGradientNormBelowTol	= impl::SolveUntilGradientNormBelowTol<true, Args...>;

// relative gradient norm below tol
template<typename ... Args>
using ConvergedWhenRelativeGradientNormBelowTol = impl::SolveUntilGradientNormBelowTol<false, Args...>;
template<typename ... Args>
using StopWhenRelativeGradientNormBelowTol	= impl::SolveUntilGradientNormBelowTol<false, Args...>;

}}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CONVERGENCE_CRITERIA_HPP_
