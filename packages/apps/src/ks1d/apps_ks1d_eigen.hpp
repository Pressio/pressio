
#ifndef ROMPPAPPS_KS1D_EIGEN_HPP_
#define ROMPPAPPS_KS1D_EIGEN_HPP_

#include "../apps_ConfigDefs.hpp"
#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include <iostream>
#include <fstream>

namespace rompp{ namespace apps{

class KS1dEigen{
  using eigVec = Eigen::VectorXd;
  using mv_t = Eigen::MatrixXd;
  using ui_t = unsigned int;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using residual_type	= eigVec;
  using jacobian_type	= Eigen::SparseMatrix
    <scalar_type, Eigen::RowMajor, int>;

  typedef Eigen::Triplet<scalar_type> Tr;

public:
  explicit KS1dEigen(eigVec params, ui_t Nnode=127)
    : mu_(params), Nnode_(Nnode){}

  KS1dEigen() = delete;
  ~KS1dEigen() = default;

public:
  void setup(){
	xR_ = mu_(1);
    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Nnode_+1);
    dxInv_ = 1.0/dx_;

    // grid
    xGrid_.resize(Nnode_);
    for (ui_t i=0; i<Nnode_; ++i)
      xGrid_(i) = dx_*(i+1);

//    // init condition
    U_.resize(Nnode_);

    // Read restart if file exists
    std::ifstream icfile("restart.inp");
    if (icfile.is_open()) {
      std::vector<scalar_type> u0;
      scalar_type u0i = 0.0;
      while (icfile >> u0i) u0.push_back(u0i);
      // Check that init_guess vector is correct length
      assert((ui_t) u0.size() == Nnode_);
      for (ui_t i=0; i < Nnode_; i++) U_(i) = u0[i];

    } else {
    	    for (ui_t i=0; i<Nnode_; ++i)
    	      U_(i) = 0.0;

    	    U_(Nnode_/2) = 1.0;
    }
    icfile.close();


    U0_ = U_;
  };

  state_type const & getInitialState() const {
    return U0_;
  };

  void residual(const state_type & u,
		residual_type & rhs,
		const scalar_type /* t */) const;

  residual_type residual(const state_type & u,
			 const scalar_type t) const{
    residual_type RR(Nnode_);
    this->residual(u, RR, t);
    return RR;
  }

  // computes: A = Jac B where B is a Eigen::MatrixXd
  void applyJacobian(const state_type & y,
		     const mv_t & B,
		     mv_t & A,
		     scalar_type t) const{
    auto JJ = jacobian(y, t);
//    std::cout << JJ << std::endl;
    // multiply
    A = JJ * B;
  }

  // computes: A = Jac B where B is a Eigen::MatrixXd
  mv_t applyJacobian(const state_type & y,
		     const mv_t & B,
		     scalar_type t) const{
    mv_t A( y.size(), B.cols() );
    applyJacobian(y, B, A, t);
    return A;
  }

  void jacobian(const state_type & u,
		jacobian_type & jac,
		const scalar_type /*t*/) const;

  jacobian_type jacobian(const state_type & u,
			 const scalar_type t) const{

    jacobian_type JJ(u.size(), u.size());
    this->jacobian(u, JJ, t);
    return JJ;
  }

private:
  eigVec mu_; // parameters
  const scalar_type xL_ = 0.0; //left side of domain
  scalar_type xR_ = 128.0; // right side of domain
  ui_t Nnode_; // # of nodes
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size
  eigVec xGrid_; // mesh points coordinates
  mutable std::vector<Tr> tripletList;
  mutable state_type U_; // state vector
  mutable state_type U0_; // initial state vector

};//end class

}} //namespace rompp::apps
#endif
