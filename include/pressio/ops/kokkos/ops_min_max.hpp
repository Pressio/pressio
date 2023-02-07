/*
//@HEADER
// ************************************************************************
//
// ops_min_max.hpp
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

#ifndef OPS_KOKKOS_OPS_MIN_MAX_HPP_
#define OPS_KOKKOS_OPS_MIN_MAX_HPP_

namespace pressio{ namespace ops{

namespace impl{

template <typename ReducerType, typename ViewType>
::pressio::mpl::enable_if_t<
  ViewType::rank == 1,
  typename ::pressio::Traits<ViewType>::scalar_type
  >
kokkos_reduce(const ViewType & view) {
  using sc_t = typename ::pressio::Traits<ViewType>::scalar_type;
  sc_t result;
  ReducerType reducer(result);
  auto f = KOKKOS_LAMBDA(size_t i, sc_t& r) {
      reducer.join(r, view(i));
    };
  Kokkos::parallel_reduce(
    "reduce_1d", view.extent(0), f, reducer);
  return result;
}

template <typename ReducerType, typename ViewType>
::pressio::mpl::enable_if_t<
  ViewType::rank == 2,
  typename ::pressio::Traits<ViewType>::scalar_type
  >
kokkos_reduce(const ViewType & view) {
  using sc_t = typename ::pressio::Traits<ViewType>::scalar_type;
  sc_t result;
  ReducerType reducer(result);
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> range_2d(
      {0, 0},
      {view.extent(0), view.extent(1)}
    );
  auto f = KOKKOS_LAMBDA(size_t i0, size_t i1, sc_t& r) {
      reducer.join(r, view(i0, i1));
    };
  Kokkos::parallel_reduce("reduce_2d", range_2d, f, reducer);
  return result;
}

}

template <typename T>
::pressio::mpl::enable_if_t<
  // TPL/container specific
    (::pressio::is_native_container_kokkos<T>::value
  || ::pressio::is_expression_acting_on_kokkos<T>::value)
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename ::pressio::Traits<T>::scalar_type
  >
max(const T & obj)
{
  assert(::pressio::ops::extent(obj, 0) > 0);

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  return impl::kokkos_reduce<Kokkos::Max<sc_t>>(impl::get_native(obj));
}

template <typename T>
::pressio::mpl::enable_if_t<
  // TPL/container specific
    (::pressio::is_native_container_kokkos<T>::value
  || ::pressio::is_expression_acting_on_kokkos<T>::value)
  // scalar compatibility
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value),
  typename ::pressio::Traits<T>::scalar_type
  >
min(const T & obj)
{
  assert(::pressio::ops::extent(obj, 0) > 0);

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  return impl::kokkos_reduce<Kokkos::Min<sc_t>>(impl::get_native(obj));
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_MIN_MAX_HPP_
