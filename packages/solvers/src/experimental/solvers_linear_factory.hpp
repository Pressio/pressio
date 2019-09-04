/*
//@HEADER
// ************************************************************************
//
// solvers_linear_factory.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP

#include <iostream>
#include <memory>
#include <type_traits>
#include <Eigen/Sparse>

#include "solvers_linear_direct.hpp"
#include "solvers_linear_iterative.hpp"
#include "solvers_linear_traits.hpp"
#include "solvers_policy_linear_iterative_trilinos.hpp"

#include "../solvers_ConfigDefs.hpp"
#include "../../../containers/src/matrix/containers_matrix_traits.hpp"


namespace pressio{ namespace solvers{



template <typename SolverT, typename MatrixT>
struct createDirectSolverTypeHelper{
  using solver_traits = linear::details::solver_traits<SolverT>;
  using wrapped_type = typename containers::details::traits<MatrixT>::wrapped_t;
  using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type>;
  using wrapped_solver_type 
    = decltype( LinearDirectSolver<concrete_solver_type, MatrixT>(std::make_shared<concrete_solver_type>()) );
};



struct LinearSolvers {

  /**
   * Create a direct linear solver for dense Eigen matrices
   *
   * @return wrapped_solver direct linear solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT,
    typename std::enable_if<
      containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createDirectSolver() 
    -> typename createDirectSolverTypeHelper<SolverT, MatrixT>::wrapped_solver_type {

    using solver_traits = linear::details::solver_traits<SolverT>;

    static_assert(solver_traits::eigen_enabled && solver_traits::direct, "Solver not available for linear systems defined by Eigen matrices");

    using wrapped_type = typename containers::details::traits<MatrixT>::wrapped_t;
    using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type>;
    //using concrete_policy_type = SolversLinearDenseEigenPolicy<concrete_solver_type, MatrixT>;

    auto solver = std::make_shared<concrete_solver_type>();
    auto wrapped_solver = LinearDirectSolver<concrete_solver_type, MatrixT>(solver);

    return wrapped_solver;
  }



template <typename SolverT, typename MatrixT>
struct createDirectSolverTypeHelper{
    using type = decltype( LinearSolvers::createDirectSolver<SolverT, MatrixT>() );
};


  /**
   * Create a direct linear solver for Eigen matrices
   *
   * @param A matrix representing the linear system to be solved
   * @return solver direct linear solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT
  >
  static auto createDirectSolver(const MatrixT& A)
    -> typename createDirectSolverTypeHelper<SolverT, MatrixT>::type{

    auto solver = LinearSolvers::createDirectSolver<SolverT, MatrixT>();
    solver.resetLinearSystem(A);
    return solver;
  }



template <typename SolverT, typename MatrixT, typename PrecT>
struct createIterativeSolverTypeHelper{
    using solver_traits = linear::details::solver_traits<SolverT>;
    using preconditioner_traits = linear::details::preconditioner_traits<PrecT>;
    using wrapped_type = typename containers::details::traits<MatrixT>::wrapped_t;
    using concrete_precon_type = typename preconditioner_traits::template eigen_preconditioner_type<wrapped_type>;
    using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type, concrete_precon_type>;
    using wrapped_solver_type 
    = decltype( LinearIterativeSolver<concrete_solver_type, MatrixT>(std::make_shared<concrete_solver_type>()) );
};

  /**
   * Create a linear iterative solver for sparse Eigen matrices
   *
   * @return wrapped_solver sparse linear iterative solver for Eigen matrices
   */
  template <
    typename SolverT,
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createIterativeSolver() 
  -> typename createIterativeSolverTypeHelper<SolverT, MatrixT, PrecT>::wrapped_solver_type {

    using solver_traits = linear::details::solver_traits<SolverT>;
    using preconditioner_traits = linear::details::preconditioner_traits<PrecT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct, "Solver not available for linear systems defined by Eigen matrices");
    static_assert(preconditioner_traits::eigen_enabled, "Preconditioner not available for linear systems defined by Eigen matrices");

    using wrapped_type = typename containers::details::traits<MatrixT>::wrapped_t;
    using concrete_precon_type = typename preconditioner_traits::template eigen_preconditioner_type<wrapped_type>;
    using concrete_solver_type = typename solver_traits::template eigen_solver_type<wrapped_type, concrete_precon_type>;
    auto solver = std::make_shared<concrete_solver_type>();
    auto wrapped_solver = LinearIterativeSolver<concrete_solver_type, MatrixT>(solver);

    return wrapped_solver;
  }



template <typename SolverT, typename MatrixT, typename PrecT>
struct createIterativeSolverTypeHelper2{
    using type = decltype( LinearSolvers::createIterativeSolver<SolverT, MatrixT>() );
};


  /**
   * Create a linear iterative solver for sparse Eigen matrices
   *
   * @param  A sparse matrix representing the linear system to be solved
   * @return solver sparse linear iterative solver for Eigen matrices
    */
  template <
    typename SolverT,
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner
  >
  static auto createIterativeSolver(const MatrixT& A)
  -> typename createIterativeSolverTypeHelper2<SolverT, MatrixT, PrecT>::type {

    auto solver = LinearSolvers::createIterativeSolver<SolverT, MatrixT>();
    solver.resetLinearSystem(A);
    return solver;
  }




//   /**
//    * Create a linear iterative solver for sparse Trilinos matrices
//    *
//    * @return An iterative linear solver for sparse Trilinos matrices
//    */
// #ifdef HAVE_TRILINOS
//   template <
//     typename SolverT,
//     typename MatrixT,
//     typename PrecT = linear::DefaultPreconditioner,
//     typename std::enable_if<
//       containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Trilinos,
//       MatrixT
//     >::type* = nullptr
//   >
//   static auto createIterativeSolver() {

//     using solver_traits = linear::details::solver_traits<SolverT>;
//     using preconditioner_traits = linear::details::preconditioner_traits<PrecT>;

//     static_assert(solver_traits::trilinos_enabled && !solver_traits::direct, "Solver not available for linear systems defined by Trilinos matrices");
//     static_assert(preconditioner_traits::trilinos_enabled, "Preconditioner not available for linear systems defined by Trilinos matrices");

//     using concrete_solver_type = AztecOO;
//     // using concrete_policy_type = SolversLinearIterativeTrilinosPolicy<concrete_solver_type, MatrixT>;

//     auto solver_flag = solver_traits::trilinos_flag;
//     auto precon_flag = preconditioner_traits::trilinos_flag;

//     auto solver = std::make_shared<concrete_solver_type>();
//     solver->SetAztecOption(AZ_solver, solver_flag);

//     switch (precon_flag) {
//       case AZ_ilu: {
//         solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
//         solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
//         break;
//       }
//       case AZ_ilut: {
//         solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
//         solver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
//         solver->SetAztecOption(AZ_overlap, 1);
//         solver->SetAztecOption(AZ_ilut_fill, 3.0);
//         break;
//       }
//       case AZ_icc: {
//         solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
//         solver->SetAztecOption(AZ_subdomain_solve, AZ_icc);
//         break;
//       }
//     }

//     auto wrapped_solver = LinearIterativeSolver<concrete_solver_type, MatrixT>(solver);

//     return wrapped_solver;
//   }

// #endif


  /**
   * Create an iterative least square linear solver for sparse Eigen matrices
   *
   * @return iterative least square linear solver for sparse Eigen matrices
   */
  template <
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createLeastSquareIterativeSolver() 
    -> decltype(LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>()) {

    return LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>();
  }


  /**
   * Create an iterative least square linear solver for sparse Eigen matrices
   *
   * @param A sparse Eigen matrix defining the least square system to be solved
   * @return iterative least square linear solver for sparse Eigen matrices
   */
  template <
    typename MatrixT,
    typename PrecT = linear::DefaultPreconditioner,
    typename std::enable_if<
      containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen,
      MatrixT
    >::type* = nullptr
  >
  static auto createLeastSquareIterativeSolver(const MatrixT& A) 
   -> decltype(LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>(A)) {

    return LinearSolvers::createIterativeSolver<typename linear::LSCG, MatrixT, PrecT>(A);
  }


};

} // end namespace solvers
} //end namespace pressio

#endif
