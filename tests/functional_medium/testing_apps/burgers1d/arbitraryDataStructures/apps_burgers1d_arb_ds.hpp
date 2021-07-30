/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_arb_ds.hpp
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

#ifndef APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_HPP_
#define APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_HPP_

#include "apps_burgers1d_arb_ds_custom_vector.hpp"
#include "apps_burgers1d_arb_ds_custom_dense_matrix.hpp"
#include <iostream>

namespace pressio{ namespace apps{

class Burgers1dArbDs
{

public:
  using int_t     = std::size_t;

  // these are detected by pressio
  using scalar_type	  = double;
  using state_type	  = arbds::Vector<scalar_type>;
  using velocity_type	  = arbds::Vector<scalar_type>;
  using jacobian_type	  = arbds::DenseMatrix<scalar_type>;

public:
  Burgers1dArbDs() = delete;

  explicit Burgers1dArbDs(int_t Ncell)
    : Ncell_(Ncell){
    setup();
  }

  int_t meshSize() const{
    return Ncell_;
  }

  state_type const & getInitialState() const{
    return U0_;
  };

public:
  velocity_type createVelocity() const{
    velocity_type f(Ncell_);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type JJ(Ncell_, Ncell_);
    return JJ;
  }

  void velocity(const state_type & u,
      const scalar_type & t,
      velocity_type & f) const{
    this->velocity_impl(u, t, f);
  }

  void jacobian(const state_type & u,
  		const scalar_type & t,
  		jacobian_type & jac) const{
    this->jacobian_impl(u, t, jac);
  }

private:
  void setup()
  {
    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    constexpr auto two = ::pressio::utils::constants<scalar_type>::two();
    constexpr auto oneHalf = one/two;
    dx_ = (xR_ - xL_)/static_cast<scalar_type>(Ncell_);
    dxInv_ = one/dx_;

    xGrid_.resize(Ncell_);
    for (int_t i=0; i<Ncell_; ++i)
      xGrid_(i) = dx_*static_cast<scalar_type>(i) + dx_*oneHalf;

    // init condition
    U_.resize(Ncell_);
    U0_.resize(Ncell_);
    for (int_t i=0; i<Ncell_; ++i){
      U_(i) = one;
      U0_(i) = one;
    }
  };

  void velocity_impl(const state_type & u,
		     const scalar_type & t,
		     velocity_type & f) const
  {
    constexpr auto one = ::pressio::utils::constants<scalar_type>::one();
    constexpr auto two = ::pressio::utils::constants<scalar_type>::two();
    constexpr auto oneHalf = one/two;

    const auto coeff = oneHalf * dxInv_;
    f(0) = coeff*(mu_[0]*mu_[0] - u(0)*u(0)) + mu_[1] * std::exp(mu_[2]*xGrid_(0));
    for (int_t i=1; i<Ncell_; ++i){
      f(i) = coeff*(u(i-1)*u(i-1) - u(i)*u(i)) + mu_[1]*std::exp(mu_[2]*xGrid_(i));
    }
  }

  void jacobian_impl(const state_type & u,
  		     const scalar_type & t,
  		     jacobian_type & jac) const
  {
    jac(0,0) = -dxInv_*u(0);
    for (int_t i=1; i<Ncell_; ++i){
      jac(i, i)   = -dxInv_ * u(i);
      jac(i, i-1) = dxInv_ * u(i-1);
    }
  }

private:
  const scalar_type xL_ = 0.0;
  const scalar_type xR_ = 100.0;
  std::array<scalar_type, 3> mu_ = {{5., 0.02, 0.02}};

  std::size_t Ncell_;
  scalar_type dx_;
  scalar_type dxInv_;
  state_type xGrid_;

  mutable state_type U_;
  state_type U0_;

};//end class

}} //namespace pressio::apps
#endif  // APPS_BURGERS1D_ARBITRARYDATASTRUCTURES_APPS_BURGERS1D_ARB_DS_HPP_
