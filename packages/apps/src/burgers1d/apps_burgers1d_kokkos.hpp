
#ifndef PRESSIOAPPS_BURGERS1D_KOKKOS_HPP_
#define PRESSIOAPPS_BURGERS1D_KOKKOS_HPP_

#include "../apps_ConfigDefs.hpp"
#include <Kokkos_Core.hpp>

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


class Burgers1dKokkos{
  using exe_space = Kokkos::DefaultExecutionSpace;
  using kll = Kokkos::LayoutLeft;

  using sc_t	= double;
  using k1dLl_d = Kokkos::View<sc_t*, kll, exe_space>;
  using k1dLl_h = k1dLl_d::HostMirror;

  using state_type_d	= k1dLl_d;
  using state_type_h	= k1dLl_h;
  using velocity_type_d	= state_type_d;
  using velocity_type_h	= state_type_h;

  using params_t = std::array<sc_t, 3>;
  using ui_t = unsigned int;

public:
  using scalar_type	= sc_t;
  using state_type	= state_type_d;
  using velocity_type	= state_type_d;

public:
  explicit Burgers1dKokkos(params_t params,
			   ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell),
      x_d_{"xd", Ncell_},
      x_h_{"xd", Ncell_},
      U_d_{"ud", Ncell_},
      U_h_{"ud", Ncell_}{}

  Burgers1dKokkos() = delete;
  ~Burgers1dKokkos() = default;

public:
  void setup(){
    dx_ = (xR_ - xL_)/static_cast<sc_t>(Ncell_);
    dxInv_ = 1.0/dx_;

    // grid
    for (ui_t i=0; i<Ncell_; ++i){
      x_h_(i) = dx_*i + dx_*0.5;
      U_h_(i) = 1.0;
    }
    Kokkos::deep_copy(x_d_, x_h_);
    Kokkos::deep_copy(U_d_, U_h_);
  };

  state_type const getInitialState() const {
    return U_d_;
  };

  void velocity(const state_type & u,
		velocity_type & rhs,
		const scalar_type /* t */) const
  {
    using func_t = VelocityFunctor<k1dLr_d, state_type, velocity_type, sc_t>;
    func_t F(mu_[0], mu_[1], mu_[2], Ncell_, dxInv_, x_d_, u, rhs);
    Kokkos::parallel_for(Ncell_, F);
  }

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const{
    velocity_type RR("RR", Ncell_);
    this->velocity(u, RR, t);
    return RR;
  }

private:
  params_t mu_; // parameters
  const sc_t xL_ = 0.0; //left side of domain
  const sc_t xR_ = 100.0; // right side of domain
  ui_t Ncell_; // # of cells
  sc_t dx_; // cell size
  sc_t dxInv_; // inv of cell size

  k1dLr_d x_d_; // mesh on device
  k1dLr_h x_h_; // mesh on host

  mutable state_type U_d_; // state on device
  mutable state_type_h U_h_; // state on host
};//end class

}} //namespace pressio::apps
#endif
