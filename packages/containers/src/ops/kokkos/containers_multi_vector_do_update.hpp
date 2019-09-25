
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include<KokkosBlas1_axpby.hpp>

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
//  overloads for computing: MV = a * MV + b * MV1
// where MV is an kokkos multivector wrapper
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & mv, const scalar_t &a,
	       const T & mv1, const scalar_t &b)
{
  KokkosBlas::axpby(b, *mv1.data(), a, *mv.data());
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & mv, const T & mv1, const scalar_t & b)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
  KokkosBlas::axpby(b, *mv1.data(), zero, *mv.data());
}

}}}//end namespace pressio::containers::ops
#endif
#endif
