
#ifndef PRESSIOAPPS_BURGERS1D_EIGEN_HPP_
#define PRESSIOAPPS_BURGERS1D_EIGEN_HPP_

#include "../apps_ConfigDefs.hpp"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <iostream>

namespace pressio{ namespace apps{

class Burgers1dEigen{
  using eigVec = Eigen::VectorXd;
  using mv_t = Eigen::MatrixXd;
  using ui_t = unsigned int;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using jacobian_type	= Eigen::SparseMatrix
    <scalar_type, Eigen::RowMajor, int>;

  typedef Eigen::Triplet<scalar_type> Tr;

public:
  explicit Burgers1dEigen(eigVec params, ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell){}

  Burgers1dEigen() = delete;
  ~Burgers1dEigen() = default;

public:
  void setup(){
    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = 1.0/dx_;

    // grid
    xGrid_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i)
      xGrid_(i) = dx_*i + dx_*0.5;

    // init condition
    U_.resize(Ncell_);
    for (ui_t i=0; i<Ncell_; ++i)
      U_(i) = 1.0;
    U0_ = U_;
  };

  state_type const & getInitialState() const {
    return U0_;
  };

  void velocity(const state_type & u,
		velocity_type & rhs,
		const scalar_type /* t */) const;

  velocity_type velocity(const state_type & u,
			 const scalar_type t) const{
    velocity_type RR(Ncell_);
    this->velocity(u, RR, t);
    return RR;
  }

  // computes: A = Jac B where B is a Eigen::MatrixXd
  void applyJacobian(const state_type & y,
		     const mv_t & B,
		     mv_t & A,
		     scalar_type t) const{
    auto JJ = jacobian(y, t);
    // std::cout << "ApplyJacobian" << std::endl;
    // std::cout << JJ << std::endl;
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
  const scalar_type xR_ = 100.0; // right side of domain
  ui_t Ncell_; // # of cells
  scalar_type dx_; // cell size
  scalar_type dxInv_; // inv of cell size
  eigVec xGrid_; // mesh points coordinates
  mutable std::vector<Tr> tripletList;
  mutable state_type U_; // state vector
  mutable state_type U0_; // initial state vector

};//end class

}} //namespace pressio::apps
#endif
