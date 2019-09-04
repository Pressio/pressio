/*
//@HEADER
// ************************************************************************
//
// solvers_converged_criterior_policy.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP
#define SOLVERS_IMPL_CONVERGE_CRITERION_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CONTAINERS_OPS"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename conv_tag>
struct IsConvergedHelper;

template <>
struct IsConvergedHelper<converged_when::completingNumMaxIters>{
  static constexpr char const * description_ = "complete max iters";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return step==maxIters;
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormCorrectionBelowTol<norm_t>>{

  static constexpr char const * description_ = "||dy|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_dy<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_r<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::relativeNormResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||R||(r) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_r/norm_r0<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::absoluteNormProjectedResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||P^T R|| < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_proj_r<tol);
  }
};

template <typename norm_t>
struct IsConvergedHelper<
  converged_when::relativeNormProjectedResidualBelowTol<norm_t>>{

  static constexpr char const * description_ = "||P^T R||(r) < tol";

  template <typename state_t, typename step_t, typename scalar_t>
  static bool evaluate(const state_t & y,
		       const state_t & dy,
		       const scalar_t & norm_dy,
		       const scalar_t & norm_r,
		       const scalar_t & norm_r0,
		       const scalar_t & norm_proj_r,
		       const scalar_t & norm_proj_r0,
		       const step_t & step,
		       const step_t & maxIters,
		       const scalar_t & tol) {
    return (norm_proj_r/norm_proj_r0<tol);
  }
};



}}}} //end namespace pressio::solvers::iterative::impl
#endif
