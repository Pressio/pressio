
#include "apps_burgers1d_eigen.hpp"

namespace pressio{ namespace apps{

void Burgers1dEigen::velocity(const state_type & u,
            velocity_type & rhs,
            const scalar_type /* t */) const{

  // std::cout << "velocity" << std::endl;
  // std::cout << rhs << std::endl;

  rhs(0) = 0.5 * dxInv_ * (mu_(0)*mu_(0) - u(0)*u(0));
  for (ui_t i=1; i<Ncell_; ++i){
    rhs(i) = 0.5 * dxInv_ * (u(i-1)*u(i-1) - u(i)*u(i));
  }
  for (ui_t i=0; i<Ncell_; ++i){
    rhs(i) += mu_(1)*exp(mu_(2)*xGrid_(i));
  }
}

void Burgers1dEigen::jacobian(const state_type & u,
            jacobian_type & jac,
            const scalar_type /*t*/)const{

  //evaluate jacobian
  if (jac.rows() == 0 || jac.cols()==0 ){
    jac.resize(u.size(), u.size());
  }
  tripletList.clear();
  tripletList.push_back( Tr( 0, 0, -dxInv_*u(0)) );
  for (ui_t i=1; i<Ncell_; ++i){
    tripletList.push_back( Tr( i, i-1, dxInv_ * u(i-1) ) );
    tripletList.push_back( Tr( i, i, -dxInv_ * u(i) ) );
  }
  jac.setFromTriplets(tripletList.begin(), tripletList.end());
}

}} //namespace pressio::apps
