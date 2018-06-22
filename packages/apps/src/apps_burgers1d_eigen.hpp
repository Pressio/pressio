
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
#include "Eigen/SparseCore"


namespace apps{

class burgers1dEigen
{
private:
  using eigVec = Eigen::VectorXd;
  using ui_t = unsigned int;
public:
  using scalar_type = double;
  using state_type = Eigen::VectorXd;//Matrix<scalar_type,-1,1>;
  using jacobian_type = Eigen::Matrix<scalar_type,-1,-1>;
  //    Eigen::SparseMatrix<scalar_type,Eigen::RowMajor,int>;

public:  
  burgers1dEigen(eigVec params, ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell){}
  ~burgers1dEigen() = default; 

  void setup(){
    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;
    // grid 
    xGrid_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i){
      xGrid_(i) = dx_*i + dx_*0.5;
    };
    // init condition
    U_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i)
      U_(i) = 1.0;
    U0_ = U_;
  };

  state_type getInitialState(){
    return U0_;
  };
  
  void residual(const state_type & u,
		state_type & rhs,
		const scalar_type /* t */)
  {
    rhs(0) = (0.5 * ( mu_(0)*mu_(0) - u(0)*u(0) ) )/dx_;
    for (ui_t i=1; i<Ncell_; ++i){
      rhs(i) = ( 0.5*(u(i-1)*u(i-1) - u(i)*u(i)) )/dx_;
    }
    for (ui_t i=0; i<Ncell_; ++i){
      rhs(i) += mu_(1)*exp(mu_(2)*xGrid_(i));
    }    
  }

  void jacobian(const state_type & u,
		jacobian_type & jac,
		const scalar_type /*t*/)
  {
    //evaluate jacobian
    if (jac.rows() == 0 || jac.cols()==0 ){
      jac.resize(u.size(), u.size());
    }
    // typedef Eigen::Triplet<double> Tr;
    // std::vector<Tr> tripletList;
    // tripletList.push_back( Tr(0,0,v_ij) );
    
    jac = jacobian_type::Zero(jac.rows(), jac.cols());
    jac(0,0) = -dxInv_ * u(0)*u(0);
    for (ui_t i=1; i<Ncell_; ++i){
      jac(i,i-1) = dxInv_ * u(i-1)*u(i-1);
      jac(i,i) = -dxInv_ * u(i-1)*u(i-1);     
    }    
  }
  
private:  
  eigVec mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain 
  const scalar_type xR_ = 100.0; // right side of domain
  ui_t Ncell_; // # of cells
  // const scalar_type t0 = 0.0; // initial time
  // const scalar_type tfinal_ = 35.0;
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size
  eigVec xGrid_; // mesh points coordinates
  state_type U_; // state vector
  state_type U0_; // initial state vector
};

}//end namespace apps
#endif 
