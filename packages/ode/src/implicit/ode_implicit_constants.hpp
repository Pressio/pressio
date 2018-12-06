
#ifndef ODE_IMPLICIT_CONSTANTS_HPP_
#define ODE_IMPLICIT_CONSTANTS_HPP_

#include "../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace coeffs{

  // beside y_n, bdf1 needs y_n-1
  static constexpr std::size_t bdf1_numAuxStates_ = 1;
  // bdf2 needs f(y_n,...)
  static constexpr std::size_t bdf1_numAuxRHS_ = 0;


  // beside y_n, bdf2 needs y_n-1, y_n-2
  static constexpr std::size_t bdf2_numAuxStates_ = 2;
  // bdf2 needs f(y_n,...)
  static constexpr std::size_t bdf2_numAuxRHS_ = 0;


  template <typename scalar_t>
  struct bdf2{
    static constexpr scalar_t c1_ = static_cast<scalar_t>(4)/3;
    static constexpr scalar_t c2_ = static_cast<scalar_t>(1)/3;
    static constexpr scalar_t c3_ = static_cast<scalar_t>(2)/3;
  };


}}}// end namespace rompp::ode::coeffs
#endif
