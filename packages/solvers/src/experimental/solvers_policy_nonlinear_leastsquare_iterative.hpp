/*
//@HEADER
// ************************************************************************
//
// solvers_policy_nonlinear_leastsquare_iterative.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_NONLINEAR_LEASTSQUARE_ITERATIVE_HPP

#include <iostream>
#include <type_traits>

#include "solvers_system_traits.hpp"
#include "solvers_linear_factory.hpp"
#include "solvers_nonlinear_traits.hpp"
#include "solvers_l2_vector_norm.hpp"
#include "solvers_meta_static_checks.hpp"
#include "../solvers_ConfigDefs.hpp"
#include "../../../CONTAINERS_OPS"
#include "../../../containers/src/meta/containers_meta_detection_idiom.hpp"


namespace pressio{
namespace solvers{

namespace todeprecate{

/**
 * Return the transpose of a matrix. To deprecate and implement in containers
 **/
template <
  typename MatrixT,
  typename ::pressio::mpl::enable_if_t<
    containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen
  >* = nullptr
>
MatrixT transpose(const MatrixT& A)
{
  MatrixT res(A.data()->transpose());
  return res;
}

} // end namespace todeprecate


/**
 * Implement the Levenberg-Marquardt algorithm for non-linear least-squares
 */
struct SolversNonLinearIterativeLeastSquareLevenbergMarquardtPolicy {

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename ::pressio::mpl::enable_if_t<
      containers::details::traits<VectorT>::is_vector &&
      solvers::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value
    >* = nullptr
  >
  static VectorT solve(
    const SystemT& sys,
    const VectorT& x0,
    containers::default_types::uint maxIterations,
    containers::default_types::uint maxNonLinearIterations,
    double tolerance,
    double nonLinearTolerance,
    double lambda
  ) {

    auto dy = sys.residual(x0);
    auto Ja = sys.jacobian(x0);
    auto JaT = todeprecate::transpose(Ja);
    auto lId = decltype(Ja)(Ja.cols(), Ja.cols());
    lId.setIdentity();
    lId.addToDiagonal(lambda-1.0);

    auto b = containers::ops::product(JaT, dy);
    auto A = containers::ops::product(JaT, Ja) + lId;

    auto solver = LinearSolvers::createIterativeSolver<SolverT, typename SystemT::matrix_type, PrecT>(A);
    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tolerance);

    auto xOld = x0;
    auto xNew = x0;

    double normO = 0.0;
    double normN = 0.0;
    containers::default_types::uint iStep = 1;

    while (iStep++ < maxNonLinearIterations) {
      xNew = xOld - solver.solve(b);
      auto dyN = sys.residual(xNew);

      normO = NormT::template compute_norm(dy);
      normN = NormT::template compute_norm(dyN);

      if (normN >= normO) {
        // Step not accepted
        lambda *= 2;
        lId.setIdentity();
        lId.addToDiagonal(lambda-1.0);

        A = containers::ops::product(JaT, Ja) + lId;
        solver.resetLinearSystem(A);
      } else {
        // Step accepted
        if ((normO-normN) < nonLinearTolerance) {break;}
        lambda *= 0.8;

        xOld = xNew;
        dy = sys.residual(xOld);
        Ja = sys.jacobian(xOld);
        JaT = todeprecate::transpose(Ja);
        lId.setIdentity();
        lId.addToDiagonal(lambda);

        b = containers::ops::product(JaT, dy);
        A = containers::ops::product(JaT, Ja) + lId;
        solver.resetLinearSystem(A);
      }
    }
    return xNew;
  }
};


/**
 * Implement the Gauss-Newton algorithm for non-linear least-squares.
 */
struct SolversNonLinearIterativeLeastSquareGaussNewtonPolicy {

  template <
    typename SolverT,
    typename PrecT,
    typename NormT,
    typename SystemT,
    typename VectorT,
    typename ::pressio::mpl::enable_if_t<
      containers::details::traits<VectorT>::is_vector &&
      solvers::meta::are_vector_compatible<
        typename details::system_traits<SystemT>::vector_type,
        VectorT
      >::value
    >* = nullptr
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
    auto JaT = todeprecate::transpose(Ja);

    auto x = x0;
    auto b = containers::ops::product(JaT, dy);
    auto A = containers::ops::product(JaT, Ja);

    auto solver = LinearSolvers::createIterativeSolver<SolverT, typename SystemT::matrix_type, PrecT>(A);
    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tolerance);

    double normN = 0.0;
    double normO = NormT::template compute_norm(dy);

    containers::default_types::uint iStep = 1;
    while (iStep++ < maxNonLinearIterations) {
      x = x - solver.solve(b);
      dy = sys.residual(x);
      normN = NormT::template compute_norm(dy);
      if (std::abs(normO - normN) < nonLinearTolerance) {break;}

      normO = normN;
      Ja = sys.jacobian(x);
      JaT = todeprecate::transpose(Ja);

      b = containers::ops::product(JaT, dy);
      A = containers::ops::product(JaT, Ja);
      solver.resetLinearSystem(A);
    }
    return x;
  }
};


} // end namespace solvers
} // end namespace pressio

#endif
