
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_TRAITS_HPP


namespace rompp{
namespace solvers{

class SolversNonLinearIterativeNewtonRaphsonPolicy; // Fwd declaration
class SolversNonLinearIterativeLeastSquareLevenbergMarquardtPolicy; // Fwd declaration

template <typename qr_algo_tag>
class SolversNonLinearIterativeLeastSquareGaussNewtonQRPolicy; // Fwd declaration


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
  using solver_type = SolversNonLinearIterativeNewtonRaphsonPolicy;
  static constexpr bool enabled = true;
};


} // end namespace details
} // end namespace nonlinear


namespace nonlinearleastsquare {

// Nonlinear least square iterative solvers types
struct LevenbergMarquardt {};
struct GaussNewtonQR {};


namespace details {

// Solvers traits
template <typename T>
struct solver_traits {
  using solver_type = void;
  static constexpr bool enabled = false;
};


template <>
struct solver_traits<LevenbergMarquardt> {
  using solver_type = SolversNonLinearIterativeLeastSquareLevenbergMarquardtPolicy;
  static constexpr bool enabled = true;
};

// template <>
// struct solver_traits<GaussNewtonQR> {
//   using solver_type = SolversNonLinearIterativeLeastSquareGaussNewtonQRPolicy;
//   static constexpr bool enabled = true;
// };


} // end namespace details
} // end namespace nonlinearleastsquare

} // end namespace solvers
} // end namespace rompp

#endif
