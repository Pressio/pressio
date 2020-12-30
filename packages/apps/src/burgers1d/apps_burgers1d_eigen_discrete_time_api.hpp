/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_eigen_discrete_time_api.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef APPS_BURGERS1D_APPS_BURGERS1D_EIGEN_DISCRETE_TIME_API_HPP_
#define APPS_BURGERS1D_APPS_BURGERS1D_EIGEN_DISCRETE_TIME_API_HPP_

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include <iostream>

namespace pressio{ namespace apps{

class Burgers1dEigenDiscreteTimeApi
{
  using eigVec = Eigen::VectorXd;
  using ui_t = unsigned int;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using discrete_time_residual_type	= eigVec;
  using jacobian_type	= Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int>;

  typedef Eigen::Triplet<scalar_type> Tr;

public:

  explicit Burgers1dEigenDiscreteTimeApi(eigVec params, ui_t Ncell=1000)
    : mu_(params), Ncell_(Ncell){
    this->setup();
  }

  Burgers1dEigenDiscreteTimeApi() = delete;
  ~Burgers1dEigenDiscreteTimeApi() = default;

public:
  state_type const & getInitialState() const {
    return U0_;
  };

  discrete_time_residual_type createDiscreteTimeResidual() const
  {
    std::cout << "calling createTimeDiscreteResidualObject" << std::endl;
    discrete_time_residual_type R(Ncell_);
    R.setConstant(0);
    return R;
  }

  Eigen::MatrixXd createApplyDiscreteTimeJacobianResult(const Eigen::MatrixXd & B) const
  {
    std::cout << "calling createApplyTimeDiscreteJacobianObject" << std::endl;
    Eigen::MatrixXd A(Ncell_, B.cols());
    A.setConstant(0);
    return A;
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    discrete_time_residual_type & R,
  			    Args && ... states) const
  {
    this->discreteTimeResidualImpl(step, time, dt,
      R,std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const Eigen::MatrixXd & B,
  				 Eigen::MatrixXd & A,
  				 Args && ... states) const
  {
    this->applyDiscreteTimeJacobianImpl(step, time, dt,
      B, A, std::forward<Args>(states)...);
  }

private:
  // case when we only have a single auxiliary state
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    const auto f =  this->velocityImpl(yn, time);
    R = yn - ynm1 - dt * f;
  }

  // case when we two auxiliary states (BDF2)
  template <typename step_t>
  void discreteTimeResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				discrete_time_residual_type & R,
				const state_type & yn,
				const state_type & ynm1,
        const state_type & ynm2) const
  {
    const auto f =  this->velocityImpl(yn, time);
    R = yn - 4./3.*ynm1 + 1./3.*ynm2 - 2./3.*dt * f;
  }


  // case when we only have a single auxiliary state
  template <typename step_t, typename state_t>
  void applyDiscreteTimeJacobianImpl(const step_t & step,
				     const scalar_type & time,
				     const scalar_type & dt,
				     const Eigen::MatrixXd & B,
				     Eigen::MatrixXd & A,
				     const state_t & yn,
				     const state_t & ynm1) const
  {
    // compute Jacobian
    auto J =  this->jacobianImpl(yn, time);

    // compute time discrete Jacobian
    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    J.coeffs() *= -dt;
    for (std::size_t i=0; i<Ncell_; ++i)
      J.coeffRef(i,i) += one;

    // compute A = J * B
    A = J * B;
  }

  // case when we have two auxiliar states (BDF2)
  template <typename step_t, typename state_t>
  void applyDiscreteTimeJacobianImpl(const step_t & step,
				     const scalar_type & time,
				     const scalar_type & dt,
				     const Eigen::MatrixXd & B,
				     Eigen::MatrixXd & A,
				     const state_t & yn,
				     const state_t & ynm1,
				     const state_t & ynm2) const
  {
    // compute Jacobian
    auto J =  this->jacobianImpl(yn, time);

    // compute time discrete Jacobian
    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    J.coeffs() *= -2./3.*dt;
    for (std::size_t i=0; i<Ncell_; ++i)
      J.coeffRef(i,i) += one;

    // compute A = J * B
    A = J * B;
  }


private:
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

  void velocityImpl(const state_type & u,
  		const scalar_type t,
		velocity_type & rhs) const
  {
    rhs(0) = 0.5 * dxInv_ * (mu_(0)*mu_(0) - u(0)*u(0));
    for (ui_t i=1; i<Ncell_; ++i){
      rhs(i) = 0.5 * dxInv_ * (u(i-1)*u(i-1) - u(i)*u(i));
    }
    for (ui_t i=0; i<Ncell_; ++i){
      rhs(i) += mu_(1)*exp(mu_(2)*xGrid_(i));
    }
  }

  velocity_type velocityImpl(const state_type & u,
  			 const scalar_type t) const{
    velocity_type RR(Ncell_);
    this->velocityImpl(u, t, RR);
    return RR;
  }

  void jacobianImpl(const state_type & u,
		const scalar_type /*t*/,
		jacobian_type & jac) const
  {
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

  jacobian_type jacobianImpl(const state_type & u,
			 const scalar_type t) const{

    jacobian_type JJ(u.size(), u.size());
    this->jacobianImpl(u, t, JJ);
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
#endif  // APPS_BURGERS1D_APPS_BURGERS1D_EIGEN_DISCRETE_TIME_API_HPP_
