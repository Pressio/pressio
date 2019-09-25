
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_SRC_OPS_KOKKOS_NORMS_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_NORMS_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "KokkosBlas1_nrm1.hpp"
#include "KokkosBlas1_nrm2.hpp"

namespace pressio{ namespace containers{ namespace ops{

template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  return KokkosBlas::nrm1(*a.data());
}


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  return KokkosBlas::nrm2(*a.data());
}

}}}//end namespace pressio::containers::ops
#endif
#endif
