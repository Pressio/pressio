
#ifndef PRESSIOAPPS_BURGERS1D_KOKKOS_FUNCTORS_HPP_
#define PRESSIOAPPS_BURGERS1D_KOKKOS_FUNCTORS_HPP_

#include "../apps_ConfigDefs.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

namespace pressio{ namespace apps{

template <
  typename x_t,
  typename y_t,
  typename rhs_t,
  typename sc_t
  >
struct VelocityFunctor{
  sc_t mu0_;
  sc_t mu1_;
  sc_t mu2_;
  int Ncell_;
  sc_t dxInv_;
  x_t x_;
  y_t y_;
  rhs_t r_;

  VelocityFunctor(sc_t mu0, sc_t mu1, sc_t mu2,
		  int Ncell, sc_t dxInv,
		  x_t x, y_t y, rhs_t r)
    : mu0_{mu0}, mu1_{mu1}, mu2_{mu2},
      Ncell_{Ncell}, dxInv_{dxInv},
      x_{x}, y_{y}, r_{r}{}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    if (i==0)
      r_(i) = 0.5 * dxInv_ * (mu0_*mu0_ - y_(i)*y_(i));

    if (i>=1 and i<Ncell_)
      r_(i) = 0.5 * dxInv_ * (y_(i-1)*y_(i-1) - y_(i)*y_(i));

    r_(i) += mu1_*std::exp(mu2_*x_(i));
  }
};


// template <
//   typename x_t,
//   typename y_t,
//   typename rhs_t,
//   typename sc_t
//   >
// struct JacobianFunctor{
//   int Ncell_;
//   sc_t dxInv_;
//   y_t y_;
//   rhs_t r_;

//   JacobianFunctor(int Ncell, sc_t dxInv,
// 		  y_t y, rhs_t r)
//     : Ncell_{Ncell}, dxInv_{dxInv},
//       y_{y}, r_{r}{}

//   KOKKOS_INLINE_FUNCTION
//   void operator() (const int & i) const {
//     // if (i==0)
//     //   r_(i) = 0.5 * dxInv_ * (mu0_*mu0_ - y_(i)*y_(i));
//     // if (i>=1 and i<Ncell_)
//     //   r_(i) = 0.5 * dxInv_ * (y_(i-1)*y_(i-1) - y_(i)*y_(i));
//     // r_(i) += mu1_*std::exp(mu2_*x_(i));
//   }
// };

}} //namespace pressio::apps
#endif
