/*
//@HEADER
// ************************************************************************
//
// rom_wls_jacobians_container_impl.hpp
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

#ifndef ROM_WLS_IMPL_ROM_WLS_JACOBIANS_CONTAINER_IMPL_HPP_
#define ROM_WLS_IMPL_ROM_WLS_JACOBIANS_CONTAINER_IMPL_HPP_

namespace pressio{ namespace rom{ namespace wls{  namespace impl{

template<typename jac_t>
class FrozenJacobiansContainer
{
  using wls_jacs_t   = std::vector<jac_t>;

public:
  FrozenJacobiansContainer(const window_size_t timeStencilSize,
			   const window_size_t numStepsInWindow,
			   const jac_t  & phi )
    :  wlsJacs_( std::min(timeStencilSize+1, numStepsInWindow)*numStepsInWindow, phi),
       jacStencilSize_(std::min(timeStencilSize+1, numStepsInWindow))
  {}

  window_size_t jacobianIndexOffset(window_size_t stepNumLocal) const{
    return stepNumLocal*jacStencilSize_;
  }

  jac_t & localJacobian(window_size_t stepNumLocal, int jacobian_index){
    return wlsJacs_[jacobianIndexOffset( stepNumLocal ) + jacStencilSize_- jacobian_index -1 ] ;
  }

  const jac_t & localJacobian(window_size_t stepNumLocal, int jacobian_index) const{
    return wlsJacs_[jacobianIndexOffset( stepNumLocal ) + jacStencilSize_- jacobian_index -1 ] ;
  }

private:
  wls_jacs_t wlsJacs_;
  window_size_t jacStencilSize_;

};


template<typename jac_t>
class NonFrozenJacobiansContainer
{
  using wls_jacs_t   = std::vector<jac_t>;

public:
  NonFrozenJacobiansContainer(const window_size_t timeStencilSize,
			      const window_size_t numStepsInWindow,
			      const jac_t  & phi)
    :  wlsJacs_( std::min(timeStencilSize+1,numStepsInWindow), phi),
       jacStencilSize_(std::min(timeStencilSize+1, numStepsInWindow))
  {}

  window_size_t jacobianIndexOffset(window_size_t stepNumLocal) const {
    return 0;
  }

  jac_t & localJacobian(window_size_t stepNumLocal, int jacobian_index) {
    return wlsJacs_[jacobianIndexOffset( stepNumLocal ) + jacStencilSize_- jacobian_index -1 ] ;
  }

  const jac_t & localJacobian(window_size_t stepNumLocal, int jacobian_index) const {
    return wlsJacs_[jacobianIndexOffset( stepNumLocal ) + jacStencilSize_- jacobian_index -1 ] ;
  }

private:
  wls_jacs_t wlsJacs_;
  window_size_t jacStencilSize_;

};

}}}} //end namespace pressio::rom::wls::impl
#endif  // ROM_WLS_IMPL_ROM_WLS_JACOBIANS_CONTAINER_IMPL_HPP_
