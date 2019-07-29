
#include "apps_burgers1d_kokkos.hpp"

#ifdef HAVE_TRILINOS
namespace pressio{ namespace apps{

void Burgers1dKokkos::setup(){
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


void Burgers1dKokkos::velocity(const state_type & u,
  const scalar_type /* t */, 
  velocity_type & rhs) const
{
  using func_t = VelocityFunctor<k1dLl_d, state_type, velocity_type, sc_t>;
  func_t F(mu_[0], mu_[1], mu_[2], Ncell_, dxInv_, x_d_, u, rhs);
  Kokkos::parallel_for(Ncell_, F);
}


typename Burgers1dKokkos::velocity_type
Burgers1dKokkos::velocity(const state_type & u,
			  const scalar_type t) const{
  velocity_type RR("RR", Ncell_);
  this->velocity(u, t, RR);
  return RR;
}


void Burgers1dKokkos::applyJacobian(const state_type & y,
				    const mvec_t & B,
				    scalar_type t,
            mvec_t & A) const{
  auto JJ = jacobian(y, t);
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();

  const char ct = 'N';
  KokkosSparse::spmv(&ct, one, JJ, B, zero, A);
}


typename Burgers1dKokkos::mvec_t
Burgers1dKokkos::applyJacobian(const state_type & y,
			       const mvec_t & B,
			       scalar_type t) const{
  mvec_t A("AA", Ncell_, B.extent(1) );
  applyJacobian(y, B, t, A);
  return A;
}

void Burgers1dKokkos::jacobian(const state_type & u,
			       const scalar_type t, 
             jacobian_type & jac) const
{
  // here the Jacobian is passed in as an argument.
  // to recompute it, use parallel for

  // state_type_h u_h = Kokkos::create_mirror_view(u);
  // jacobian_type JJ2("JJ2", numRows, numCols, numEnt, val, ptr, ind);
  // jac = JJ2;
  //Kokkos::deep_copy(jac.values, JJ2.values);

  using exe_space = typename Burgers1dKokkos::execution_space;

  using policy_type = Kokkos::RangePolicy<exe_space, int>;

  using func_t = JacobianFunctor<state_type, jacobian_type, sc_t>;
  func_t F(Ncell_, dxInv_, u, jac);
  Kokkos::parallel_for( policy_type(0, Ncell_), F);
  //Kokkos::parallel_for(Ncell_, F);
}


typename Burgers1dKokkos::jacobian_type
Burgers1dKokkos::jacobian(const state_type & u,
			  const scalar_type t) const
{
  const int numRows = Ncell_;
  const int numCols = Ncell_;
  const int numEnt = (Ncell_-1)*2 + 1;

  // here we need to create a Jacobian, so use the following.
  // The data is filled on the host, then copied to device

  state_type_h u_h = Kokkos::create_mirror_view(u);

  typename jacobian_type::row_map_type::non_const_type ptr ("ptr", numRows+1);
  {
    typename jacobian_type::row_map_type::HostMirror ptr_h = Kokkos::create_mirror_view (ptr);
    ptr_h[0] = 0;
    for (int lclRow = 0; lclRow < numRows; ++lclRow) {
      auto nElPerRow = (lclRow==0) ? 1 : 2;
      ptr_h[lclRow+1] = ptr_h[lclRow] + nElPerRow; // 1 entry in row 0, 2 in all others
    }
    Kokkos::deep_copy (ptr, ptr_h);
  }

  typename jacobian_type::index_type::non_const_type ind ("ind", numEnt);
  {
    typename jacobian_type::index_type::HostMirror ind_h = Kokkos::create_mirror_view (ind);
    ind_h[0] = 0;
    int k = 0;
    for (int lclRow = 1; lclRow < numRows; ++lclRow) {
      ind_h[++k] = lclRow-1;
      ind_h[++k] = lclRow;
    }
    Kokkos::deep_copy (ind, ind_h);
  }

  typename jacobian_type::values_type val ("val", numEnt);
  {
    typename jacobian_type::values_type::HostMirror val_h = Kokkos::create_mirror_view (val);
    val_h[0] = -dxInv_*u_h(0);
    int k = 0;
    for (int lclRow = 1; lclRow < numRows; ++lclRow) {
      val_h[++k] = dxInv_*u_h[lclRow-1];
      val_h[++k] = -dxInv_*u_h[lclRow];
    }
    Kokkos::deep_copy (val, val_h);
  }

  jacobian_type JJ("JJ", numRows, numCols, numEnt, val, ptr, ind);
  this->jacobian(u, t, JJ);
  return JJ;
}

}} //namespace pressio::apps
#endif
