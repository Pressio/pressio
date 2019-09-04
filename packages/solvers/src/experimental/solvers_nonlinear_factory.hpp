/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear_factory.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP

#include "../solvers_ConfigDefs.hpp"
#include "solvers_nonlinear_iterative.hpp"
#include "solvers_nonlinear_leastsquare_iterative.hpp"
#include "solvers_nonlinear_traits.hpp"
#include "solvers_policy_nonlinear_iterative.hpp"
#include "solvers_policy_nonlinear_leastsquare_iterative.hpp"
#include "solvers_l2_vector_norm.hpp"
#include "solvers_linear_traits.hpp"
#include "../../../containers/src/meta/containers_meta_detection_idiom.hpp"
#include "../../../QR_BASIC"


namespace pressio {
namespace solvers {

struct NonLinearSolvers {

  /**
   * Raise an exception while trying to create an invalid nonlinear solver.
   */
  template <
    typename NSolverT,
		typename LSolverT,
    typename std::enable_if<
      !nonlinear::details::solver_traits<NSolverT>::enabled,
      void
    >::type* = nullptr
  >
  static void createIterativeSolver() {
  	std::cerr << "Error: the nonlinear solver selected \
is not available or its name was mispelt" << std::endl;
  	assert(0);
  }
  //--------------------------------------------------------------


  template <typename NSolverT, typename LSolverT>
  struct createIterativeSolverTypeHelper{
    using solver_traits = linear::details::solver_traits<LSolverT>;
    using policy_type = typename nonlinear::details::solver_traits<NSolverT>::solver_type;
    using ret_type = decltype( NonLinearIterativeSolver<policy_type, LSolverT>() );

  };


  /**
   * Create a nonlinear solver.
   */
  template <
    typename NSolverT,
		typename LSolverT,
    typename std::enable_if<
      nonlinear::details::solver_traits<NSolverT>::enabled,
      void
    >::type* = nullptr
  >
  static auto createIterativeSolver()
  -> typename createIterativeSolverTypeHelper<NSolverT, LSolverT>::ret_type {

    using solver_traits = linear::details::solver_traits<LSolverT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct,
		  "Error: either the linear solver is a direct one \
or is not available for linear systems defined by Eigen matrices");

    using policy_type =
      typename nonlinear::details::solver_traits<NSolverT>::solver_type;
    return NonLinearIterativeSolver<policy_type, LSolverT>();
  }
  //--------------------------------------------------------------


  template <typename NSolverT, typename LSolverT>
  struct createNonLinIterativeLSSolverTypeHelper{
    using solver_traits = linear::details::solver_traits<LSolverT>;
    using policy_type = typename nonlinearleastsquare::details::solver_traits<NSolverT>::solver_type;
    using ret_type = NonLinearLeastSquareIterativeSolver<policy_type, LSolverT>;

  };

  /**
   * Create a nonlinear least square iterative solver
   */
  template <
    typename NSolverT,
    typename LSolverT,
    typename ::pressio::mpl::enable_if_t<
      nonlinearleastsquare::details::solver_traits<NSolverT>::enabled
    >* = nullptr
  >
  static auto createNonLinearIterativeLeastSquareSolver()
  -> typename createNonLinIterativeLSSolverTypeHelper<NSolverT, LSolverT>::ret_type {

    using solver_traits = linear::details::solver_traits<LSolverT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct,
		  "Error: either the linear solver is a direct one \
or is not available for linear systems defined by Eigen matrices");

    using policy_type = typename nonlinearleastsquare::details::solver_traits<NSolverT>::solver_type;
    return NonLinearLeastSquareIterativeSolver<policy_type, LSolverT>();
  }
  //--------------------------------------------------------------

  // /**
  //  * Create a nonlinear least square iterative solver
  //  * using QR for each inner linear least square solve
  //  */
  // template <
  //   typename NSolverT,
  //   typename qr_algo_tag = ::pressio::qr::Hacked,
  //   typename ::pressio::mpl::enable_if_t<
  //     std::is_same<NSolverT, nonlinearleastsquare::GaussNewtonQR>::value
  //   >* = nullptr
  // >
  // static SolversNonLinearIterativeLeastSquareGaussNewtonQRPolicy<qr_algo_tag> 
  // createNonLinearIterativeLeastSquareQRBasedSolver(){
  //   return SolversNonLinearIterativeLeastSquareGaussNewtonQRPolicy<qr_algo_tag>();
  // }
  // //--------------------------------------------------------------


};

} // end namespace solvers
} // end namespace pressio
#endif
