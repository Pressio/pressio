/*
//@HEADER
// ************************************************************************
//
// solvers_admissible_qr_solver_for_gn_qr.hpp
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

#ifndef SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_QR_SOLVER_FOR_GN_QR_HPP_
#define SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_QR_SOLVER_FOR_GN_QR_HPP_

namespace pressio{ namespace nonlinearsolvers{

#ifdef PRESSIO_ENABLE_CXX20
template <class T, class StateType, class MatrixType, class RType>
concept QRSolverForGnQr =
  /* compound requirements */
  requires(T & A,
	   const MatrixType & matrix,
	   const RType & operand,
	   StateType & b,
	   const StateType & d)
  {
    { A.computeThin(matrix) } -> std::same_as<void>;
    { A.applyQTranspose(operand, b) } -> std::same_as<void>;
    { A.applyRTranspose(d, b) } -> std::same_as<void>;
    { A.solve(d, b) } -> std::same_as<void>;
  };
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::nonlinearsolvers

#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace nonlinearsolvers{

template <class T, class StateType, class MatrixType, class RType, class enable = void>
struct QRSolverForGnQr
  : std::false_type{};

template <class T, class StateType, class MatrixType, class RType>
struct QRSolverForGnQr<
  T, StateType, MatrixType, RType,
  ::pressio::mpl::enable_if_t<
  	// must have computeThin
	std::is_void<
	  decltype
	  (
	   std::declval<T>().computeThin
	   (std::declval<MatrixType const &>())
	   )
	  >::value and
  	// must have applyQTranspose
	std::is_void<
	  decltype
	  (
	   std::declval<T const>().applyQTranspose
	   (
	    std::declval<RType const &>(),
	    std::declval<StateType &>()
	    )
	   )
	  >::value and
  	// must have applyRTranspose
	std::is_void<
	  decltype
	  (
	   std::declval<T const>().applyRTranspose
	   (
	    std::declval<StateType const &>(),
	    std::declval<StateType &>()
	    )
	   )
	  >::value and
  	// must have solve
	std::is_void<
	  decltype
	  (
	   std::declval<T const>().solve
	   (
	    std::declval<StateType const &>(),
	    std::declval<StateType &>()
	    )
	   )
	  >::value
	>
  > : std::true_type{};

}} // end namespace pressio::nonlinearsolvers
#endif

#endif  // SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_QR_SOLVER_FOR_GN_QR_HPP_
