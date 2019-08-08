
#ifndef PRESSIOAPPS_BURGERS1D_KOKKOS_FUNCTORS_HPP_
#define PRESSIOAPPS_BURGERS1D_KOKKOS_FUNCTORS_HPP_

#include "../apps_ConfigDefs.hpp"

#ifdef HAVE_KOKKOS
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

namespace pressio{ namespace apps{

template <
  typename x_t,
  typename u_t,
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
  u_t u_;
  rhs_t r_;

  VelocityFunctor(sc_t mu0, sc_t mu1, sc_t mu2,
		  int Ncell, sc_t dxInv,
		  x_t x, u_t y, rhs_t r)
    : mu0_{mu0}, mu1_{mu1}, mu2_{mu2},
      Ncell_{Ncell}, dxInv_{dxInv},
      x_{x}, u_{y}, r_{r}{}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    if (i==0)
      r_(i) = 0.5 * dxInv_ * (mu0_*mu0_ - u_(i)*u_(i));

    if (i>=1 and i<Ncell_)
      r_(i) = 0.5 * dxInv_ * (u_(i-1)*u_(i-1) - u_(i)*u_(i));

    r_(i) += mu1_*std::exp(mu2_*x_(i));
  }
};


template <
  typename u_t,
  typename crs_mat_t,
  typename sc_t
  >
struct JacobianFunctor{
  int Ncell_;
  sc_t dxInv_;
  u_t u_;
  crs_mat_t J_;

  JacobianFunctor(int Ncell,
		  sc_t dxInv,
		  u_t y,
		  crs_mat_t Jin)
    : Ncell_{Ncell},
      dxInv_{dxInv},
      u_{y},
      J_{Jin}{}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & iRow) const
  {
    const int ncol = (iRow==0) ? 1 : 2;
    int cols[2];
    sc_t vals[2];

    if (iRow==0){
      cols[0] = 0;
      vals[0] = -dxInv_*u_(iRow);
    }
    else{
      cols[0] = iRow-1;
      vals[0] = dxInv_*u_(iRow-1);

      cols[1] = iRow;
      vals[1] = -dxInv_*u_(iRow);
    }

    J_.replaceValues(iRow, cols, ncol, vals, false, true);
  }
};

}} //namespace pressio::apps
#endif
#endif
