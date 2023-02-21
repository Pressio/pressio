/*
//@HEADER
// ************************************************************************
//
// solvers_admissible_state.hpp
//                          Pressio
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

#ifndef PRESSIO_CONCEPTS_VECTOR_SPACE_ELEMENT_HPP_
#define PRESSIO_CONCEPTS_VECTOR_SPACE_ELEMENT_HPP_

namespace pressio{

#ifdef PRESSIO_ENABLE_CXX20

template <class T>
concept VectorSpaceElement =
  requires(){
    // the trait specialization AND scalar_type have to exists
    typename ::pressio::Traits<T>::scalar_type;
  }
  /*
  axiom ScalarTypeRepresentsTheFieldScalar(){
    // scalar_type represents the scalar type of the field F on which the space is defined
  } &&
  axiom AdditiveClosure(const T& v, const T& w){
    // v + w is defined and belongs to same vector space
  } &&
  axiom ScalarMultiplicationClosure(const T& v,
				    typename pressio::Traits<T>::scalar_type a){
    // av is defined and belongs to same space
  } &&
  axiom Associativity(const T& u, const T& v, const T& w){
    // u + (v + w) = (u + v) + w
  } &&
  axiom Commutativity(const T& u, const T& v){
    // u + v = v + u
  } &&
  axiom IdentityElement(const T& v){
    // there exits a zero-element "0 \in V" such that: v + 0 = v
  } &&
  axiom InverseElement(const T& v){
    // there exists an element -v such that: v + (−v) = 0 (the zero-element)
  } &&
  axiom ScalarCompatibility(const T& v,
			    typename pressio::Traits<T>::scalar_type a,
			    typename pressio::Traits<T>::scalar_type b){
    // a(bv) = (ab)v
  } &&
  axiom ScalarMultiplicationIdentity(const T& v){
    // 1v = v where 1 is the multiplicative identity in F
  } &&
  axiom DitributivityWrtVectorAddition(const T& u,
				       const T& v,
				       typename pressio::Traits<T>::scalar_type a,
				       typename pressio::Traits<T>::scalar_type b){
    // a(u + v) = au + av
  } &&
  axiom DitributivityWrtFieldAddition(const T& v,
				      typename pressio::Traits<T>::scalar_type a,
				      typename pressio::Traits<T>::scalar_type b){
    // (a + b)v = av + bv
  }*/
  ;



// template<typename T, typename enable = void>
// struct VectorSpaceElement : std::false_type{};

// template<typename T>
// struct VectorSpaceElement<
//   T,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::has_scalar_typedef<
//       ::pressio::Traits<T>
//       >::value
//     >
//   > : std::true_type{};


// template<typename T>
// using RealVectorSpaceElement =  VectorSpaceElement<T>;

} // namespace pressio
#endif
