
#ifndef SOLVERS_EXPERIMENTAL_TRAITS_HPP_
#define SOLVERS_EXPERIMENTAL_TRAITS_HPP_

#include <Eigen/Core>


namespace solvers {


// Forward declarations
class NonlinearNewtonRaphsonSolver;


namespace linear {

// Linear iterative solvers types
struct CG {};
struct Gmres {};
struct Bicgstab {};

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

  template <typename T>
  using eigen_solver_type = Eigen::ConjugateGradient<T>;

  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct solver_traits<Bicgstab> {

  template <typename T>
  using eigen_solver_type = Eigen::BiCGSTAB<T>;

  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};


// Preconditioners traits 
template <typename T>
struct preconditioner_traits {
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = false;
};

template<>
struct preconditioner_traits<DefaultPreconditioner> {
  static constexpr bool eigen_enabled = true;
  static constexpr bool trilinos_enabled = true;
};

template <>
struct preconditioner_traits<Jacobi> {
  static constexpr bool eigen_enabled = false;
  static constexpr bool trilinos_enabled = true;
};

} // end namespace details
} // end namespace linear


namespace nonlinear {

// Nonlinear iterative solvers types
struct NewtonRaphson {};


namespace details {

// Solvers traits
template <typename T>
struct solver_traits {
  
  using solver_type = void;

  static constexpr bool solver_exists = false; 
};


template <>
struct solver_traits<NewtonRaphson> {

  using solver_type = NonlinearNewtonRaphsonSolver;
  
  static constexpr bool solver_exists = true;
};


} // end namespace details
} // end namespace nonlinear
} // end namespace solvers

#endif
