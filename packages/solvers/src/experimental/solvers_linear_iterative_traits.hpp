#ifndef SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_ITERATIVE_TRAITS_HPP

#include <Eigen/Core>
#include "AztecOO.h"

#include "matrix/core_matrix_traits.hpp"


namespace solvers {
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
  using eigen_solver_type = Eigen::ConjugateGradient<MatrixT, Eigen::Lower, PrecT>;

  static constexpr int trilinos_flag = AZ_cg;

  static constexpr bool dense_only = false;
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<Gmres> {

  static constexpr int trilinos_flag = AZ_gmres;

  static constexpr bool dense_only = false;
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<Bicgstab> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::BiCGSTAB<MatrixT, PrecT>;

  static constexpr int trilinos_flag = AZ_bicgstab;

  static constexpr bool dense_only = false;
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<ColPivHouseholderQR> {

  template <
    typename MatrixT
  >
  using eigen_solver_type = Eigen::ColPivHouseholderQR<MatrixT>;

  static constexpr bool dense_only = true;
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<CompleteOrthogonalDecomposition> {

  template <
    typename MatrixT
  >
  using eigen_solver_type = Eigen::CompleteOrthogonalDecomposition<MatrixT>;

  static constexpr bool dense_only = true;
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<LSCG> {

  template <
    typename MatrixT,
    typename PrecT = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>
  >
  using eigen_solver_type = Eigen::LeastSquaresConjugateGradient<MatrixT, PrecT>;

  static constexpr bool dense_only = false;
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = false;
};


// Preconditioners traits
template <typename T>
struct preconditioner_traits {
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = false;
};

template<>
struct preconditioner_traits<DefaultPreconditioner> {

  template <typename MatrixT>
  using eigen_preconditioner_type = Eigen::DiagonalPreconditioner<typename MatrixT::Scalar>;

  static constexpr int trilinos_flag = INT_MIN;

  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct preconditioner_traits<Jacobi> {
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = false;
};

} // end namespace details
} // end namespace linear
} // end namespace solvers

#endif
