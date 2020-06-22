/*
//@HEADER
// ************************************************************************
//
// solvers_norm_dispatcher.hpp
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

#ifndef SOLVERS_IMPL_NORM_DISPATCHER_HPP
#define SOLVERS_IMPL_NORM_DISPATCHER_HPP

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ud_ops_t>
struct NormDispatcher
{
private:
  const ud_ops_t * udOps_ = nullptr;

public:
  NormDispatcher() = default;
  NormDispatcher(const ud_ops_t * udOps) : udOps_{udOps}{}

  template < typename vec_t, typename scalar_t >
  mpl::enable_if_t<
    !::pressio::containers::meta::is_vector_wrapper_arbitrary<vec_t>::value
    >
  evaluate(const vec_t & vecIn,
	   scalar_t & result,
	   const ::pressio::solvers::Norm & normType) const
  {
    if (normType == ::pressio::solvers::Norm::L1){
      result = ::pressio::ops::norm1(vecIn);
    }
    else if (normType == ::pressio::solvers::Norm::L2){
      result = ::pressio::ops::norm2(vecIn);
    }
    else
      throw std::runtime_error("Invalid norm type, cannot conmpute norm");
  }


  template < typename vec_t, typename scalar_t, typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_arbitrary<vec_t>::value and
    !std::is_void<_ud_ops_t>::value
    >
  evaluate(const vec_t & vecIn,
	   scalar_t & result,
	   const ::pressio::solvers::Norm & normType) const
  {
    if (normType == ::pressio::solvers::Norm::L1){
      result = udOps_->norm1( *vecIn.data() );
    }
    else if (normType == ::pressio::solvers::Norm::L2){
      result = udOps_->norm2( *vecIn.data() );
    }
    else
      throw std::runtime_error("Invalid norm type, cannot conmpute norm");
  }
};

}}}} //end namespace pressio::solvers::iterative::impl
#endif
