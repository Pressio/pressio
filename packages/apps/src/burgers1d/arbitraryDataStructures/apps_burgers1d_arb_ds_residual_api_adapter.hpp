/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_arb_ds_residual_api_adapter.hpp
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

#ifndef PRESSIOAPPS_BURGERS1D_ARB_DS_RESIDUAL_API_ADAPTER_HPP_
#define PRESSIOAPPS_BURGERS1D_ARB_DS_RESIDUAL_API_ADAPTER_HPP_

#include "../apps_ConfigDefs.hpp"
#include "apps_burgers1d_arb_ds.hpp"

namespace pressio{ namespace apps{

class Burgers1dArbDsResidualApiAdapter
{
  using int_t		  = int32_t;
  using sc_t		  = double;

  using scalar_type	  = sc_t;
  using state_type	  = arbds::Vector<scalar_type>;
  using residual_type	  = arbds::Vector<scalar_type>;
  using jacobian_type	  = arbds::DenseMatrix<scalar_type>;
  using dense_matrix_type = arbds::DenseMatrix<scalar_type>;

public:
  Burgers1dArbDsResidualApiAdapter() = delete;

  Burgers1dArbDsResidualApiAdapter(const Burgers1dArbDs & appObj)
    : appObj_{appObj}
  {}

public:
  template <typename step_t, typename ... Args>
  void timeDiscreteResidual(const step_t & step,
  			    const scalar_type & time,
  			    const scalar_type & dt,
  			    residual_type & R,
  			    Args && ... states) const
  {
    //timeDiscreteResidualImpl(step, time, dt, R, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void applyTimeDiscreteJacobian(const step_t & step,
  				 const scalar_type & time,
  				 const scalar_type & dt,
  				 const dense_matrix_type & B,
  				 int id,
  				 dense_matrix_type & A,
  				 Args && ... states) const
  {
    //applyTimeDiscreteJacobianImpl(step, time, dt, B, id, A, std::forward<Args>(states)...);
  }

  residual_type createTimeDiscreteResidualObject(const state_type & stateIn) const
  {
    std::cout << "calling createTimeDiscreteResidualObject" << std::endl;
    residual_type R(Ncell_);
    //R.setConstant(0);
    return R;
  }

  dense_matrix_type createApplyTimeDiscreteJacobianObject(const state_type & stateIn,
							  const dense_matrix_type & B) const
  {
    std::cout << "calling createApplyTimeDiscreteJacobianObject" << std::endl;
    //dense_matrix_type A(Ncell_, B.cols());
    //A.setConstant(0);
    return A;
  }

private:
  // case when we only have a single auxiliary state
  template <typename step_t>
  void timeDiscreteResidualImpl(const step_t & step,
				const scalar_type & time,
				const scalar_type & dt,
				residual_type & R,
				const state_type & yn,
				const state_type & ynm1) const
  {
    // const auto f =  this->velocity(yn, time);
    // R = yn - ynm1 - dt * f;
  }

  // case when we only have a single auxiliary state
  template <typename step_t, typename state_t>
  void applyTimeDiscreteJacobianImpl(const step_t & step,
				     const scalar_type & time,
				     const scalar_type & dt,
				     const dense_matrix_type & B,
				     int id,
				     dense_matrix_type & A,
				     const state_t & yn,
				     const state_t & ynm1) const
  {
    // compute Jacobian
    auto J =  this->jacobian(yn, time);

    // // compute time discrete Jacobian
    // constexpr auto one = ::pressio::utils::constants::one<scalar_type>();
    // J.coeffs() *= -dt;
    // for (std::size_t i=0; i<Ncell_; ++i)
    //   J.coeffRef(i,i) += one;

    // // compute A = J * B
    // A = J * B;
  }

private:
  const Burgers1dArbDs & appObj_;
};//end class

}} //namespace pressio::apps
#endif
