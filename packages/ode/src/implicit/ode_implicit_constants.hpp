
#ifndef ODE_IMPLICIT_CONSTANTS_HPP_
#define ODE_IMPLICIT_CONSTANTS_HPP_

#include "../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{ namespace coeffs{

  template <typename scalar_t>
  struct bdf2{
    static constexpr scalar_t c1 = static_cast<scalar_t>(4.)/3.;
    static constexpr scalar_t c2 = static_cast<scalar_t>(1.)/3.;
    static constexpr scalar_t c3 = static_cast<scalar_t>(2.)/3.;
  };

}}}}// end namespace rompp::ode::impl::coeffs
#endif
