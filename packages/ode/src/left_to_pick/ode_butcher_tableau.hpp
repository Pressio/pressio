
#ifndef ODE_BUTCHER_TABLEAU_HPP_
#define ODE_BUTCHER_TABLEAU_HPP_

#include "ode_ConfigDefs.hpp"
#include <Eigen/Dense>

namespace ode{
namespace impl{
  
template <typename scalar, int nr, int nc>
class butcherTableau{
public:
using a_t = Eigen::Array<scalar, nr, nc>;
using c_t = Eigen::Array<scalar, nr, 1>;
using b_t = Eigen::Array<scalar, 1, nc>;
public:
  butcherTableau()
    : a_(a_t::Zero()),
      b_(b_t::Zero()),
      c_(c_t::Zero())
  {}
public:
  scalar c(int i) const{
    return c_[i];
  }
  scalar a(int i, int j) const{
    return a_(i,j);
  }
  scalar b(int i) const{
    return b_[i];
  }
  
protected:
  a_t a_;
  b_t b_;
  c_t c_;
};//end class

  
}//end namespace impl
}//end namespace ode  
#endif

