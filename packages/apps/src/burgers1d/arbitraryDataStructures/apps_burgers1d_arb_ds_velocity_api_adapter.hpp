/*
//@HEADER
// ************************************************************************
//
// apps_burgers1d_arb_ds_velocity_api_adapter.hpp
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

#ifndef PRESSIOAPPS_BURGERS1D_ARB_DS_VELOCITY_API_ADAPTER_HPP_
#define PRESSIOAPPS_BURGERS1D_ARB_DS_VELOCITY_API_ADAPTER_HPP_

#include "../../apps_ConfigDefs.hpp"
#include "apps_burgers1d_arb_ds.hpp"

namespace pressio{ namespace apps{

class Burgers1dArbDsVelocityApiAdapter
{
public:
  using int_t		  = typename Burgers1dArbDs::int_t;
  using scalar_type	  = typename Burgers1dArbDs::scalar_type;
  using state_type	  = typename Burgers1dArbDs::state_type;
  using velocity_type    = typename Burgers1dArbDs::velocity_type;
  using residual_type	  = typename Burgers1dArbDs::state_type;
  using dense_matrix_type = typename Burgers1dArbDs::dense_matrix_type;
  using jacobian_type	  = typename Burgers1dArbDs::dense_matrix_type;

  Burgers1dArbDsVelocityApiAdapter() = delete;

  Burgers1dArbDsVelocityApiAdapter(const Burgers1dArbDs & appObj)
    : appObj_{appObj},
      f_{appObj.meshSize()},
      JJ_{appObj.meshSize(), appObj.meshSize()}{}

public:
  void velocity(const state_type & u,
  		const scalar_type & t,
  		velocity_type & f) const{
    appObj_.velocity(u, t, f);
  }

  velocity_type velocity(const state_type & u,
			 const scalar_type & t) const{
    appObj_.velocity(u, t, f_);
    return f_;
  }

  // computes: A = Jac B
  void applyJacobian(const state_type & y,
		     const dense_matrix_type & B,
		     scalar_type t,
		     dense_matrix_type & A) const{
    appObj_.jacobian(y, t, JJ_);

    for (int_t i=0; i<A.extent(0); ++i){
      for (int_t j=0; j<A.extent(1); ++j){
	A(i,j) = ::pressio::utils::constants::zero<scalar_type>();
	for (int_t k=0; k<JJ_.extent(1); ++k){
	  A(i,j) += JJ_(i,k) * B(k,j);
	}
      }
    }
  }

  // computes: A = Jac B
  dense_matrix_type applyJacobian(const state_type & y,
				  const dense_matrix_type & B,
				  scalar_type t) const{
    dense_matrix_type A( y.extent(0), B.extent(1) );
    this->applyJacobian(y, B, t, A);
    return A;
  }

private:
  const Burgers1dArbDs & appObj_;
  mutable state_type f_;
  mutable jacobian_type JJ_;

};//end class

}} //namespace pressio::apps
#endif
