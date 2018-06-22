
#ifndef APPS_BURGERS1D_STDVECTOR_HPP_
#define APPS_BURGERS1D_STDVECTOR_HPP_

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


namespace apps{


class burgers1dStdVector
{
private:
  using ui_t = unsigned int;

public:
  using state_type = std::vector<double>;

public:  
  burgers1dStdVector(std::vector<double> params) : mu_(params)
  {}

  void setup(){
    dx_ = (xR_ - xL_)/static_cast<double>(Ncell_);
    xGrid_.resize(Ncell_,0.0);
    for (ui_t i=0; i<Ncell_; ++i){
      xGrid_[i] = dx_*i + dx_*0.5;
    };
    U_.resize(Ncell_, 1.0);
  };

  state_type getInitialState(){
    return U_;
  };
  
  void operator() ( const state_type & u,
		    state_type & R,
		    const double /* t */ )
  {
    R[0] = (0.5 * ( mu_[0]*mu_[0] - u[0]*u[0] ) )/dx_;
    for (ui_t i=1; i<Ncell_; ++i){
      R[i] = ( 0.5*(u[i-1]*u[i-1] - u[i]*u[i]) )/dx_;
    }
    for (ui_t i=0; i<Ncell_; ++i){
      R[i] += mu_[1]*exp(mu_[2]*xGrid_[i]);
    }    
  }

private:  
  const double xL_ = 0.0;
  const double xR_ = 100.0;
  const ui_t Ncell_ = 10;
  std::vector<double> mu_;
  const double t0 = 0.0;
  const double tfinal_ = 35.0;
  double dx_;
  std::vector<double> xGrid_;
  state_type U_;
};

  

}//end namespace apps

#endif 
