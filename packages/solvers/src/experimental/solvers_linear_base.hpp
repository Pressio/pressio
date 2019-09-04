/*
//@HEADER
// ************************************************************************
//
// solvers_linear_base.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP

#include <memory>
#include <type_traits>

#include "solvers_meta_static_checks.hpp"


namespace pressio{
namespace solvers{


/**
 * Base class for linear solver implemented through CRTP.
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear solver class.
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class LinearSolvers.
 */
template<
  typename SolverT,
  typename MatrixT,
  typename Derived
>
class LinearSolverBase {

public:

  /**
   * Initialize a new linear solver.
   *
   * @param A matrix representing the linear system to be solved.
   */
  template <
  typename CompatibleMatrixT,
  typename std::enable_if<
    solvers::meta::are_matrix_compatible<
      MatrixT,
      CompatibleMatrixT
      >::value, MatrixT*
    >::type = nullptr
  >
  void resetLinearSystem(const CompatibleMatrixT& A) {
    solver_->resetLinearSystem(A);
  }


  /**
   * Solve the linear system
   *
   * @param b RHS vector
   * @return solution vector
   */
  template <
    typename VectorLT,
    typename std::enable_if<
      solvers::meta::are_vector_matrix_compatible<
	VectorLT,
	MatrixT
        >::value, MatrixT*
      >::type = nullptr
    >
  VectorLT solve(const VectorLT& b){
    return this->underlying()._solve(b);
  }


  /**
   * Specify and solve the linear system
   *
   * @param A matrix representing the linear system to be solved
   * @param b RHS vector
   * @return solution vector
   */
  template <
    typename CompatibleMatrixT,
    typename VectorRT,
    typename std::enable_if<
      solvers::meta::are_vector_matrix_compatible<
	VectorRT,
	CompatibleMatrixT
        >::value,
      CompatibleMatrixT*
      >::type = nullptr
    >
  auto solve(const CompatibleMatrixT& A, const VectorRT& b)
    -> decltype(this->solve(b)){
    this->resetLinearSystem(A);
    return this->solve(b);
  }


  /**
   * @brief  Solve the linear system
   *
   * @param  b is the RHS vector
   * @param  x is the solution vector
   * @return void
   */
  template <
    typename VectorLT,
    typename VectorRT,
    typename std::enable_if<
      solvers::meta::are_vector_compatible<
	VectorLT,
	VectorRT
        >::value,
      VectorLT*
      >::type = nullptr
    >
  void solve(const VectorLT& b, VectorRT& x) {
    x = VectorRT(*this->solve(b).data());
  }


  /**
   * @brief  Specify and solve the linear system
   *
   * @param  A is the system matrix
   * @param  b is the RHS vector
   * @param  x is the solution vector
   * @return void
   */
  template <
    typename CompatibleMatrixT,
    typename VectorLT,
    typename VectorRT
    >
  void solve(const CompatibleMatrixT& A, const VectorLT& b, VectorRT& x) {
    this->resetLinearSystem(A);
    this->solve(b, x);
  }


protected:

  LinearSolverBase() : solver_(nullptr) {};


  LinearSolverBase(std::shared_ptr<SolverT> solver) : solver_(solver) {}


  LinearSolverBase(LinearSolverBase&& other) : solver_(std::move(other.solver_)) {}


  LinearSolverBase(const LinearSolverBase&) = delete;


  virtual ~LinearSolverBase() = default;


  std::shared_ptr<SolverT> getSolver() {
    return solver_;
  }


private:

  Derived& underlying() {
    return static_cast<Derived&>(*this);
  }


  Derived const& underlying() const {
    return static_cast<Derived const&>(*this);
  }


private:

  std::shared_ptr<SolverT> solver_;
};

} //end namespace solvers

}//end namespace pressio
#endif
