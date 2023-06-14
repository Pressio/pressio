/*
//@HEADER
// ************************************************************************
//
// ops_scale.hpp
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

#ifndef OPS_EIGEN_OPS_SCALE_HPP_
#define OPS_EIGEN_OPS_SCALE_HPP_

namespace pressio{ namespace ops{

/* constrained via is_convertible because the impl is using
   Eigen native operations which are based on expressions and require
   coefficients to be convertible to scalar types of the vector/matrix operand */
template <typename T, class ScalarType>
::pressio::mpl::enable_if_t<
  // rank-1 update common constraints
    (::pressio::Traits<T>::rank == 1
  || ::pressio::Traits<T>::rank == 2)
  // TPL/container specific
  && (::pressio::is_native_container_eigen<T>::value
  || ::pressio::is_expression_acting_on_eigen<T>::value)
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<ScalarType, typename ::pressio::Traits<T>::scalar_type>::value
  >
scale(T & o, const ScalarType & value)
{
  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  constexpr sc_t zero = ::pressio::utils::Constants<sc_t>::zero();
  sc_t value_(value);
  if (value_ == zero) {
    ::pressio::ops::set_zero(o);
  } else {
    auto && on = impl::get_native(o);
    on *= value_;
  }
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_SCALE_HPP_
