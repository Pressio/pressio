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

#include "../../../algebra/src/matrix/algebra_matrix_traits.hpp"


namespace rompp{
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
}//end namespace rompp
#endif
