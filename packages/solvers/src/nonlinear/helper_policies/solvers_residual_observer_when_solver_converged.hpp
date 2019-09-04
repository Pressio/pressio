/*
//@HEADER
// ************************************************************************
//
// solvers_residual_observer_when_solver_converged.hpp
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

#ifndef SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP
#define SOLVERS_IMPL_RESIDUAL_OBSERVER_WHEN_SOLVER_CONVERGED_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../meta/solvers_is_legitimate_residual_observer_when_converged.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template<typename observer_t,
	 typename residual_t,
	 typename = void>
struct ResidualObserverWhenSolverConverged{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid) {
    // no op
  }
};


template<typename observer_t,
	 typename residual_t>
struct ResidualObserverWhenSolverConverged<
  observer_t, residual_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::solvers::meta::is_legitimate_residual_observer_when_solver_converged<
      observer_t, residual_t
     >::value
    >
  >{

  static void evaluate(const observer_t * obs,
		       const residual_t & resid) {
    obs->observeResidualWhenSolverConverged(resid);
  }
};


}}}} //end namespace pressio::solvers::iterative::impl
#endif
