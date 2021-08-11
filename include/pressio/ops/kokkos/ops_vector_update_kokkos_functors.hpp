/*
//@HEADER
// ************************************************************************
//
// ops_vector_update_kokkos_functors.hpp
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

#ifndef OPS_KOKKOS_OPS_VECTOR_UPDATE_KOKKOS_FUNCTORS_HPP_
#define OPS_KOKKOS_OPS_VECTOR_UPDATE_KOKKOS_FUNCTORS_HPP_

namespace pressio{ namespace ops{ namespace impl{

template <class T1, class T2, class T3, class sc_t>
struct DoUpdateTwoTermsFunctor {
  T1 v_ = {};
  T2 v1_ = {};
  T3 v2_ = {};
  sc_t a_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t b_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t c_ = ::pressio::utils::Constants<sc_t>::zero();

  DoUpdateTwoTermsFunctor(T1 v, T2 v1, T3 v2,
			  sc_t a, sc_t b, sc_t c)
    : v_{v}, v1_{v1}, v2_{v2},
      a_{a}, b_{b}, c_{c}{}

  DoUpdateTwoTermsFunctor(T1 v, T2 v1, T3 v2,
			  sc_t b, sc_t c)
    : v_{v}, v1_{v1}, v2_{v2},
      b_{b}, c_{c}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i);
  }
};

template <class T1, class T2, class T3, class T4, class sc_t>
struct DoUpdateThreeTermsFunctor {
  T1 v_ = {};
  T2 v1_ = {};
  T3 v2_ = {};
  T4 v3_ = {};
  sc_t a_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t b_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t c_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t d_ = ::pressio::utils::Constants<sc_t>::zero();

  DoUpdateThreeTermsFunctor(T1 v, T2 v1, T3 v2, T4 v3,
			    sc_t a, sc_t b, sc_t c, sc_t d)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3},
      a_{a}, b_{b}, c_{c}, d_{d}{}

  DoUpdateThreeTermsFunctor(T1 v, T2 v1, T3 v2, T4 v3,
			    sc_t b, sc_t c, sc_t d)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3},
      b_{b}, c_{c}, d_{d}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i) + d_*v3_(i);
  }
};


template <class T1, class T2, class T3, class T4, class T5, class sc_t>
struct DoUpdateFourTermsFunctor {
  T1 v_ = {};
  T2 v1_ = {};
  T3 v2_ = {};
  T4 v3_ = {};
  T5 v4_ = {};
  sc_t a_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t b_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t c_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t d_ = ::pressio::utils::Constants<sc_t>::zero();
  sc_t e_ = ::pressio::utils::Constants<sc_t>::zero();

  DoUpdateFourTermsFunctor(T1 v, T2 v1, T3 v2, T4 v3, T5 v4,
			   sc_t a, sc_t b, sc_t c, sc_t d, sc_t e)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3}, v4_{v4},
      a_{a}, b_{b}, c_{c}, d_{d}, e_{e}{}

  DoUpdateFourTermsFunctor(T1 v, T2 v1, T3 v2, T4 v3, T5 v4,
			   sc_t b, sc_t c, sc_t d, sc_t e)
    : v_{v}, v1_{v1}, v2_{v2}, v3_{v3}, v4_{v4},
      b_{b}, c_{c}, d_{d}, e_{e}{}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
    v_(i) = a_*v_(i) + b_*v1_(i) + c_*v2_(i) + d_*v3_(i) + e_*v4_(i);
  }
};


}}}//end namespace pressio::ops::impl
#endif  // OPS_KOKKOS_OPS_VECTOR_UPDATE_KOKKOS_FUNCTORS_HPP_
