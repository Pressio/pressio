
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_TRAITS_HPP


namespace solvers {

class NonlinearNewtonRaphsonSolver; // Fwd declaration


namespace nonlinear {

// Nonlinear iterative solvers types
struct NewtonRaphson {};


namespace details {

// Solvers traits
template <typename T>
struct solver_traits {
  using solver_type = void;
  static constexpr bool enabled = false; 
};


template <>
struct solver_traits<NewtonRaphson> {
  using solver_type = NonlinearNewtonRaphsonSolver;
  static constexpr bool enabled = true;
};


} // end namespace details
} // end namespace nonlinear
} // end namespace solvers

#endif