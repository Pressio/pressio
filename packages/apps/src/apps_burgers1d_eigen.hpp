

#ifndef APPS_BURGERS1D_EIGEN_HPP_
#define APPS_BURGERS1D_EIGEN_HPP_

#include "apps_ConfigDefs.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>
#include "Eigen/Dense"


namespace apps{


class burgers1dEigen
{
private:
  using eigVec = Eigen::VectorXd;
  using ui_t = unsigned int;

public:
  using state_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;

public:  
  burgers1dEigen(eigVec params, ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell)
  {}

  void setup(){
    dx_ = (xR_ - xL_)/static_cast<double>(Ncell_);
    xGrid_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i){
      xGrid_(i) = dx_*i + dx_*0.5;
    };
    U_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i)
      U_(i) = 1.0;
    U0_ = U_;
  };

  state_type copyInitialState(){ return U0_; };
  
  void operator() ( const state_type & u,
		    state_type & rhs,
		    const double /* t */ )
  {
    rhs(0) = (0.5 * ( mu_(0)*mu_(0) - u(0)*u(0) ) )/dx_;
    for (ui_t i=1; i<Ncell_; ++i){
      rhs(i) = ( 0.5*(u(i-1)*u(i-1) - u(i)*u(i)) )/dx_;
    }
    for (ui_t i=0; i<Ncell_; ++i){
      rhs(i) += mu_(1)*exp(mu_(2)*xGrid_(i));
    }    
  }

  void operator() ( const state_type & u,
		    state_type & rhs,
		    jacobian_type & jac,
		    const double t )
  {
    (*this)(u,rhs,t);

    double dxInv = 1.0/dx_;
    //evaluate jacobian
    jac.resize(u.size(), u.size());
    jac = jacobian_type::Zero(jac.rows(), jac.cols());
    
    jac(0,0) = -dxInv * u(0)*u(0);
    for (ui_t i=1; i<Ncell_; ++i){
      jac(i,i-1) = dxInv * u(i-1)*u(i-1);
      jac(i,i) = -dxInv * u(i-1)*u(i-1);     
    }    
  }

  
private:  
  const double xL_ = 0.0;
  const double xR_ = 100.0;
  ui_t Ncell_;
  eigVec mu_;
  const double t0 = 0.0;
  const double tfinal_ = 35.0;
  double dx_;
  eigVec xGrid_;
  state_type U_;
  state_type U0_;
};

}//end namespace apps
#endif 
