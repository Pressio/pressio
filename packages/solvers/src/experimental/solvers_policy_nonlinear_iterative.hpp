/*
//@HEADER
// ************************************************************************
//
// solvers_policy_nonlinear_iterative.hpp
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
#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_ITERATIVE_HPP

#include <iostream>
#include <type_traits>

#include "solvers_system_traits.hpp"
#include "../solvers_ConfigDefs.hpp"
#include "solvers_linear_factory.hpp"
#include "solvers_meta_static_checks.hpp"


namespace pressio {
namespace solvers {


struct SolversNonLinearIterativeNewtonRaphsonPolicy {

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename std::enable_if<
      !containers::details::traits<VectorT>::is_vector ||
      !solvers::meta::are_vector_compatible<
          typename details::system_traits<SystemT>::vector_type,
          VectorT
        >::value,
      int
    >::type* = nullptr
  >
  static VectorT solve(
    const SystemT& system,
    const VectorT& x0,
    containers::default_types::uint maxIterations,
    containers::default_types::uint maxNonLinearIterations,
    double tolerance,
    double nonLinearTolerance
  ) {

    std::cerr << "Error: the type of the RHS vector \
is not compatible with the provided nonlinear system" << std::endl;
    assert(0);

  	return x0;
  }
  //--------------------------------------------------------------

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename std::enable_if<
      containers::details::traits<VectorT>::is_vector &&
      solvers::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value,
      int
    >::type* = nullptr
  >
  static VectorT solve(
    const SystemT& sys,
    const VectorT& x0,
    containers::default_types::uint maxIterations,
    containers::default_types::uint maxNonLinearIterations,
    double tolerance,
    double nonLinearTolerance
  ) {

    auto dy = sys.residual(x0);
    auto Ja = sys.jacobian(x0);

    auto solver = LinearSolvers::createIterativeSolver<
      SolverT, typename SystemT::matrix_type, PrecT>(Ja);
    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tolerance);

    containers::default_types::uint iStep = 1;
    VectorT xOld = x0;
    VectorT xNew(x0 - solver.solve(dy));

    while (iStep++ < maxNonLinearIterations &&
	   NormT::template compute_norm_difference(xOld, xNew)
	   > nonLinearTolerance) {
        xOld = xNew;
        sys.residual(xNew, dy);
        sys.jacobian(xNew, Ja);

        solver.resetLinearSystem(Ja);
        xNew -= solver.solve(dy);
    }
    return xNew;
  }
  //--------------------------------------------------------------

};

} // end namespace solvers
} //end namespace pressio

#endif
