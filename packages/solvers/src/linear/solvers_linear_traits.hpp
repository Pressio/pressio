#ifndef SOLVERS_EXPERIMENTAL_LINEAR_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_TRAITS_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_tags.hpp"
#include "../../../core/src/matrix/core_matrix_traits.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Core>

// #ifdef HAVE_ARMADILLO
//   #include "solvers_linear_wrapper_armadillo.hpp"
// #endif
// #ifdef HAVE_TRILINOS
//   #include "AztecOO.h"
// #endif


namespace rompp{ namespace solvers{ namespace linear { namespace details {

// Solvers traits
template <typename T>
struct traits {
  static constexpr bool direct        = false;
  static constexpr bool iterative     = false;
  static constexpr bool eigen_enabled = false;
};


template <>
struct traits<CG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type             = Eigen::ConjugateGradient<MatrixT, Eigen::Lower, PrecT>;

  static constexpr bool direct        = false;
  static constexpr bool eigen_enabled = true;
};


template <>
struct traits<Bicgstab> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type             = Eigen::BiCGSTAB<MatrixT, PrecT>;

  static constexpr bool direct        = false;
  static constexpr bool eigen_enabled = true;
};


template <>
struct traits<LSCG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::LeastSquaresConjugateGradient<MatrixT, PrecT>;

  static constexpr bool direct = false;
  static constexpr bool eigen_enabled = true;
};



// template <>
// struct traits<ColPivHouseholderQR> {

//   template <typename MatrixT>
//   using eigen_solver_type = Eigen::ColPivHouseholderQR<MatrixT>;

//   static constexpr bool direct = true;
//   static constexpr bool eigen_enabled = true;
// };


// template <>
// struct traits<Gmres> {

//   static constexpr bool direct = false;
//   static constexpr bool eigen_enabled = false;
// };


// template <>
// struct traits<CompleteOrthogonalDecomposition> {

//   template <
//     typename MatrixT
//   >
//   using eigen_solver_type = SolversLinearDirectWrapperEigen<Eigen::CompleteOrthogonalDecomposition<MatrixT>>;

//   static constexpr bool direct = true;
//   static constexpr bool eigen_enabled = true;
// };




}}}}//end namespace rompp::solvers::linear::details
#endif
