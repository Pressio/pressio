/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_binary_expression_templates.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_VECTOR_VECTOR_SHAREDMEM_BINARY_EXPRESSION_TEMPLATES_HPP_
#define CONTAINERS_VECTOR_VECTOR_SHAREDMEM_BINARY_EXPRESSION_TEMPLATES_HPP_

#include "../containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace exprtemplates{


template <typename der_t>
class SharedMemVecExpressionBase {
  SharedMemVecExpressionBase() = default;
  ~SharedMemVecExpressionBase() = default;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<der_t>::type;

public:
  static constexpr bool is_shared_mem = true;
  static constexpr bool is_vector_expression = true;
};
///-----------------------------------------------------

template <typename T,
	  typename enable = void>
struct is_sharedmem_vector_expression : std::false_type{};

template <typename T>
struct is_sharedmem_vector_expression<T,
      ::pressio::mpl::enable_if_t<
       ::pressio::mpl::publicly_inherits_from<
	T,SharedMemVecExpressionBase<T>
	 >::value
	>> : std::true_type{};


//----------------------------------------------------
// default
template <typename OP_t, typename T1,
	  typename T2, typename value_t,
	  typename ord_t, typename enable = void>
class SharedMemVectorBinaryExp
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>>{

  const T1 & a_ = {};
  const T2 & b_ = {};
  OP_t op_ = {};
  using this_t = SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>;
  friend SharedMemVecExpressionBase<this_t>;

 public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(T1 const& a, T2 const& b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_(i));}

  ord_t size() const {
    return a_.size();}
};


//----------------------------------------------------
// T1: whatever, T2: vector
template <typename OP_t, typename T1, typename T2,
	  typename value_t, typename ord_t>
class SharedMemVectorBinaryExp<
         OP_t, T1, T2, value_t, ord_t,
	 ::pressio::mpl::enable_if_t<
	   !std::is_scalar<T1>::value &&
	   containers::meta::is_vector_wrapper<T2>::value &&
	   containers::details::traits<T2>::is_shared_mem>
  >
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>>{

  const T1 & a_ = {};
  const T2 & b_ = {};
  OP_t op_ = {};
  using this_t = SharedMemVectorBinaryExp<OP_t,T1,T2,value_t,ord_t>;
  friend SharedMemVecExpressionBase<this_t>;

public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(T1 const& a, T2 const& b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_(i));}

  ord_t size() const {
    return b_.size();}
};


//-----------------------------------------------------
// T1: not scalar, T2: scalar
template <typename OP_t, typename T1,
	    typename value_t, typename ord_t>
class SharedMemVectorBinaryExp<
         OP_t, T1, value_t, value_t, ord_t,
	 ::pressio::mpl::enable_if_t<
	   !std::is_scalar<T1>::value &&
	   std::is_scalar<value_t>::value
	   > >
  : public SharedMemVecExpressionBase<
  SharedMemVectorBinaryExp<OP_t,T1,value_t,value_t,ord_t>>{

  const T1 & a_ = {};
  value_t b_ = {};
  OP_t op_ = {};
  using th_t = SharedMemVectorBinaryExp<OP_t,T1,value_t,value_t,ord_t>;
  friend SharedMemVecExpressionBase<th_t>;

public:
  using sc_type = value_t;
  using ord_type = ord_t;

  SharedMemVectorBinaryExp(const T1 & a, value_t b)
    : a_(a), b_(b){}
  ~SharedMemVectorBinaryExp() = default;

  value_t operator()(size_t i) const {
    return op_(a_(i), b_);}

  ord_t size() const{
    return a_.size();}
};


}}}//end namespace pressio::containers::exprtemplates
#endif
