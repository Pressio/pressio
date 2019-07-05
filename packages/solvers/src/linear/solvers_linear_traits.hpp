#ifndef SOLVERS_EXPERIMENTAL_LINEAR_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_TRAITS_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_tags.hpp"
#include "../../../containers/src/matrix/containers_matrix_traits.hpp"
#include "../../../containers/src/matrix/containers_matrix_meta.hpp"
#include "../../../containers/src/multi_vector/containers_multi_vector_meta.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Core>

// #ifdef HAVE_ARMADILLO
//   #include "solvers_linear_wrapper_armadillo.hpp"
// #endif
// #ifdef HAVE_TRILINOS
//   #include "AztecOO.h"
// #endif


namespace pressio{ namespace solvers{ namespace linear { namespace details {

// Solvers traits
template <typename T>
struct traits {
  static constexpr bool direct        = false;
  static constexpr bool iterative     = false;
  static constexpr bool eigen_enabled = false;
  static constexpr bool kokkos_enabled = false;
};


template <>
struct traits<::pressio::solvers::linear::iterative::CG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::ConjugateGradient<MatrixT, Eigen::Lower, PrecT>;

  static constexpr bool direct        = false;
  static constexpr bool eigen_enabled = true;
};


template <>
struct traits<::pressio::solvers::linear::iterative::Bicgstab> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::BiCGSTAB<MatrixT, PrecT>;

  static constexpr bool direct        = false;
  static constexpr bool eigen_enabled = true;
};


template <>
struct traits<::pressio::solvers::linear::iterative::LSCG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::LeastSquaresConjugateGradient<MatrixT, PrecT>;

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = true;
};


template <>
struct traits<::pressio::solvers::linear::direct::ColPivHouseholderQR> {

  template <typename MatrixT>
  using eigen_solver_type = Eigen::ColPivHouseholderQR<MatrixT>;

  static constexpr bool direct = true;
  static constexpr bool eigen_enabled = true;
};

#ifdef HAVE_TRILINOS
template <>
struct traits<::pressio::solvers::linear::direct::getrs> {

  static constexpr bool direct = true;
  // for now, disable eigen, enable it later
  static constexpr bool eigen_enabled = false;
  static constexpr bool kokkos_enabled = true;
};
#endif


// template <>
// struct traits<CompleteOrthogonalDecomposition> {

//   template <typename MatrixT>
//   using eigen_solver_type = Eigen::CompleteOrthogonalDecomposition<MatrixT>;

//   static constexpr bool direct = true;
//   static constexpr bool eigen_enabled = true;
// };


}}}}//end namespace pressio::solvers::linear::details
#endif
