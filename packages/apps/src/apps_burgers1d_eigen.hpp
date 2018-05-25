

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

public:  
  burgers1dEigen(eigVec params) : mu_(params)
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
  };

  state_type getInitialState(){
    return U_;
  };
  
  void operator() ( const state_type & u,
		    state_type & R,
		    const double /* t */ )
  {
    R(0) = (0.5 * ( mu_(0)*mu_(0) - u(0)*u(0) ) )/dx_;
    for (ui_t i=1; i<Ncell_; ++i){
      R(i) = ( 0.5*(u(i-1)*u(i-1) - u(i)*u(i)) )/dx_;
    }
    for (ui_t i=0; i<Ncell_; ++i){
      R(i) += mu_(1)*exp(mu_(2)*xGrid_(i));
    }    
  }

private:  
  const double xL_ = 0.0;
  const double xR_ = 100.0;
  const ui_t Ncell_ = 100;
  eigVec mu_;
  const double t0 = 0.0;
  const double tfinal_ = 35.0;
  double dx_;
  eigVec xGrid_;

  state_type U_;
};

  
}//end namespace apps
#endif 
