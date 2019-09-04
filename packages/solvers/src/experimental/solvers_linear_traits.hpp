/*
//@HEADER
// ************************************************************************
//
// solvers_linear_traits.hpp
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
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_TRAITS_HPP

#include <Eigen/Core>
#include "../solvers_ConfigDefs.hpp"
#include "solvers_linear_wrapper_eigen.hpp"

#ifdef HAVE_ARMADILLO
  #include "solvers_linear_wrapper_armadillo.hpp"
#endif

#ifdef HAVE_TRILINOS
  #include "AztecOO.h"
#endif

#include "../../../containers/src/matrix/containers_matrix_traits.hpp"


namespace pressio{
namespace solvers{
namespace linear {

// Linear dense solvers types
struct ColPivHouseholderQR {};
struct CompleteOrthogonalDecomposition {};

// Linear iterative solvers types
struct CG {};
struct Gmres {};
struct Bicgstab {};
struct LSCG {};

// Preconditioner types
struct Jacobi {};
struct DefaultPreconditioner {};


namespace details {

// Solvers traits
template <typename T>
struct solver_traits {
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = false;
};

template <>
struct solver_traits<CG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = SolversLinearIterativeWrapperEigen<Eigen::ConjugateGradient<MatrixT, Eigen::Lower, PrecT>>;

#ifdef HAVE_TRILINOS
  static constexpr int trilinos_flag = AZ_cg;
#endif

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = true;
#endif
};


template <>
struct solver_traits<Gmres> {

#ifdef HAVE_TRILINOS
  static constexpr int trilinos_flag = AZ_gmres;
#endif

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = false;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = true;
#endif
};


template <>
struct solver_traits<Bicgstab> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = SolversLinearIterativeWrapperEigen<Eigen::BiCGSTAB<MatrixT, PrecT>>;

#ifdef HAVE_TRILINOS
  static constexpr int trilinos_flag = AZ_bicgstab;
#endif

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = true;
#endif
};

template <>
struct solver_traits<ColPivHouseholderQR> {

  template <
    typename MatrixT
  >
  using eigen_solver_type = SolversLinearDirectWrapperEigen<Eigen::ColPivHouseholderQR<MatrixT>>;

  static constexpr bool direct = true;
  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = false;
#endif
};

template <>
struct solver_traits<CompleteOrthogonalDecomposition> {

  template <
    typename MatrixT
  >
  using eigen_solver_type = SolversLinearDirectWrapperEigen<Eigen::CompleteOrthogonalDecomposition<MatrixT>>;

  static constexpr bool direct = true;
  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = false;
#endif
};

template <>
struct solver_traits<LSCG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = SolversLinearIterativeWrapperEigen<Eigen::LeastSquaresConjugateGradient<MatrixT, PrecT>>;

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = false;
#endif
};


// Preconditioners traits
template <typename T>
struct preconditioner_traits {
  static constexpr bool eigen_enabled = false;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = false;
#endif
};

template<>
struct preconditioner_traits<DefaultPreconditioner> {

  template <typename MatrixT>
  using eigen_preconditioner_type = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>;

  static constexpr int trilinos_flag = INT_MIN;

  static constexpr bool eigen_enabled = true;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = true;
#endif
};

template <>
struct preconditioner_traits<Jacobi> {
  static constexpr bool eigen_enabled = false;
#ifdef HAVE_TRILINOS
  static constexpr bool trilinos_enabled = false;
#endif
};

} // end namespace details
} // end namespace linear
} // end namespace solvers
}//end namespace pressio
#endif
